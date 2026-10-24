#!/usr/bin/env bash
# make_cohort_fixture.sh — regenerate the cohort integration fixture deterministically.
#
# Provenance record of the committed tests/data/cohort/ fixtures. It builds them through a
# FASTQ -> BAM -> VCF upstream that is not part of noncodingdrivers (the pipeline starts
# from called variants), so it needs the tools below outside the pipeline.
#
# Produces a 5-tumour cohort on real GRCh38 chr17 sequence with a known planted
# recurrent driver (chr17:7,675,028 C>T at a CpG site, element REG_DRIVER) plus
# 60 unique passenger SNVs, then runs it through the shared upstream to per-sample
# VCFs. This is the reproducible source of the `test_cohort` Nextflow profile and
# of the harmonized fixtures committed under tests/data/cohort/.
#
# Requires: samtools, bcftools, bwa-mem2, fastp, wgsim, gatk4 (BioContainers or a
# conda env), plus a GRCh38 chr17 reference bundle (GRCh38_chr17.fa.gz and its bwa-mem2 index).
#
# Usage: benchmarks/make_cohort_fixture.sh <refbundle_dir> <out_dir>
set -euo pipefail

REFDIR="${1:?path to GRCh38 chr17 reference bundle}"
OUT="${2:?output dir}"
mkdir -p "$OUT"; cd "$OUT"

REF="$REFDIR/GRCh38_chr17.fa.gz"
INDEX="$(find -L "$REFDIR" -name '*.bwt.2bit.64' | sed 's/.bwt.2bit.64//')"

REGION_START=7660000
REGION_END=7690000
DRIVER_POS=7675028          # CpG C>T, high functional impact + recurrent
N_SAMPLES=5
SEED=7                       # controls passenger placement; keep fixed for reproducibility

# 1. Build mutant references (recurrent driver in all tumours + per-sample passengers)
#    and a pure-WT matched normal. Deterministic under $SEED.
python3 - "$REFDIR" "$REGION_START" "$REGION_END" "$DRIVER_POS" "$N_SAMPLES" "$SEED" <<'PY'
import sys, gzip, os, random
refdir, s, e, drv, n, seed = sys.argv[1:]
s, e, drv, n, seed = int(s), int(e), int(drv), int(n), int(seed)
seq=[]
with gzip.open(f"{refdir}/GRCh38_chr17.fa.gz","rt") as fh:
    fh.readline()
    for line in fh: seq.append(line.strip())
chr17="".join(seq)
region=chr17[s-1:e]; loc=drv-s
assert region[loc].upper()=='C' and region[loc+1].upper()=='G', "driver site is not a CpG"
random.seed(seed); COMP=set("ACGT")
def wfa(p,nm,x):
    open(p,"w").write(f">{nm}\n"+"\n".join(x[i:i+60] for i in range(0,len(x),60))+"\n")
gt=open("cohort_ground_truth.tsv","w"); gt.write("chrom\tpos\tref\talt\tclass\tn_samples\n")
gt.write(f"chr17\t{drv}\tC\tT\trecurrent_driver_CpG\t{n}\n")
for k in range(1,n+1):
    x=list(region); x[loc]='T'; picks=0
    while picks<12:
        p=random.randint(200,len(region)-200)
        if p==loc: continue
        b=region[p].upper()
        if b not in COMP: continue
        alt=random.choice([y for y in "ACGT" if y!=b]); x[p]=alt
        gt.write(f"chr17\t{s+p}\t{b}\t{alt}\tpassenger_T{k}\t1\n"); picks+=1
    wfa(f"T{k}_mut.fa","region","".join(x))
wfa("normal.fa","region",region); gt.close()
print(f"planted driver chr17:{drv} C>T (CpG) in all {n} tumours + {n*12} passengers")
PY

# 2. Matched normal: align + markdup once
wgsim -N 10000 -1 150 -2 150 -e 0.001 -r 0 -R 0 -S 99 normal.fa NORM_R1.fq NORM_R2.fq >/dev/null 2>&1
bwa-mem2 mem -t 4 -R "@RG\tID:NORMAL\tSM:NORMAL\tPL:ILLUMINA" "$INDEX" NORM_R1.fq NORM_R2.fq 2>/dev/null \
  | samtools sort -@4 -o NORMAL.sorted.bam -
gatk MarkDuplicates -I NORMAL.sorted.bam -O NORMAL.md.bam -M NORMAL.metrics.txt 2>/dev/null
samtools index NORMAL.md.bam

# 3. Per-tumour: reads (65% WT + 35% mutant -> ~35% VAF) -> fastp -> bwa-mem2 -> markdup -> Mutect2 -> norm
: > cohort_samplesheet.csv; echo "sample,fastq_1,fastq_2" >> cohort_samplesheet.csv
for k in $(seq 1 $N_SAMPLES); do
  wgsim -N 6500 -1 150 -2 150 -e 0.001 -r 0 -R 0 -S $((k*10+1)) normal.fa   w1.fq w2.fq >/dev/null 2>&1
  wgsim -N 3500 -1 150 -2 150 -e 0.001 -r 0 -R 0 -S $((k*10+2)) T${k}_mut.fa m1.fq m2.fq >/dev/null 2>&1
  cat w1.fq m1.fq > T${k}_R1.fq; cat w2.fq m2.fq > T${k}_R2.fq
  fastp -i T${k}_R1.fq -I T${k}_R2.fq -o T${k}.t1.fq.gz -O T${k}.t2.fq.gz --json T${k}.fastp.json --thread 4 2>/dev/null
  bwa-mem2 mem -t 4 -R "@RG\tID:T${k}\tSM:T${k}\tPL:ILLUMINA" "$INDEX" T${k}.t1.fq.gz T${k}.t2.fq.gz 2>/dev/null \
    | samtools sort -@4 -o T${k}.sorted.bam -
  gatk MarkDuplicates -I T${k}.sorted.bam -O T${k}.md.bam -M T${k}.metrics.txt 2>/dev/null
  samtools index T${k}.md.bam
  gatk Mutect2 -R "$REF" -I T${k}.md.bam -I NORMAL.md.bam -normal NORMAL \
    -L chr17:${REGION_START}-${REGION_END} -O T${k}.raw.vcf.gz 2>/dev/null
  gatk FilterMutectCalls -R "$REF" -V T${k}.raw.vcf.gz -O T${k}.filt.vcf.gz 2>/dev/null
  bcftools norm -m -any -f "$REF" T${k}.filt.vcf.gz 2>/dev/null \
    | bcftools view -f PASS,. -Oz -o T${k}.norm.vcf.gz
  bcftools index -t -f T${k}.norm.vcf.gz
  echo "T${k},T${k}.t1.fq.gz,T${k}.t2.fq.gz" >> cohort_samplesheet.csv
done
rm -f w1.fq w2.fq m1.fq m2.fq T*_R1.fq T*_R2.fq NORM_R*.fq

echo "cohort fixture ready in $OUT"
echo "  per-sample normalized VCFs: T{1..$N_SAMPLES}.norm.vcf.gz"
echo "  pool these -> OncodriveFML mutation file / ActiveDriverWGS mutation table,"
echo "  run both callers, harmonize -> the committed tests/data/cohort/ fixtures."
