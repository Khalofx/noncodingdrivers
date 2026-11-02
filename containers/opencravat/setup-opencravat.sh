#!/usr/bin/env bash
# Shared by the Dockerfile and opencravat.def, run once at build time after the conda and pip
# locks are installed. 1) installs the OpenCRAVAT modules of modules.freeze.json (exact versions)
# from the store with a build-time configuration whose metrics address is already the closed local
# port; 2) installs the run-time configuration cravat-system.yml (every remote address local) and
# leaves the metrics folder empty; 3) checks the result.
set -euo pipefail
build=/opt/opencravat/build
pkg=$(python -c "import importlib.util, os; print(os.path.dirname(importlib.util.find_spec('cravat').origin))")
mkdir -p "$pkg/conf" /opt/opencravat/modules /opt/opencravat/metrics /opt/opencravat/workspace /opt/opencravat/jobs /opt/opencravat/logs

# 1. build-time configuration: the store is reachable, the metrics address is not
sed -e 's|^store_url: .*|store_url: https://store.opencravat.org|' "$build/cravat-system.yml" > "$pkg/conf/cravat-system.yml"
cat "$pkg/cravat-system.template.yml" | sed -n '/^base_modules:/,$p' >> "$pkg/conf/cravat-system.yml"
oc module install-freeze "$build/modules.freeze.json" --yes
oc module freeze > "$build/installed.freeze.json"
python - "$build/modules.freeze.json" "$build/installed.freeze.json" <<'PY'
import json, sys
want = {(m["name"], m["version"]) for m in json.load(open(sys.argv[1]))}
have = {(m["name"], m["version"]) for m in json.load(open(sys.argv[2]))}
assert want == have, f"installed modules differ from modules.freeze.json: missing {sorted(want - have)}, extra {sorted(have - want)}"
print(f"{len(have)} OpenCRAVAT modules installed at the frozen versions")
PY

# UCSC liftOver chain files are licensed for non-commercial use only: none stays in the image.
# open-cravat bundles two (hg19ToHg38, hg18ToHg38) for `-l hg19|hg18` input, which the pipeline
# never uses (the OPENCRAVAT module fixes -l hg38); the hg19 annotator, which ships a third, is
# not installed (encode_tfbs, its only user, is not in modules.freeze.json).
rm -f "$pkg"/liftover/*.over.chain
chains=$(find / -xdev \( -name '*.over.chain' -o -name '*.over.chain.gz' \) 2>/dev/null || true)
if [ -n "$chains" ]; then echo "UCSC chain files left in the image: $chains" >&2; exit 1; fi

# 2. run-time configuration: every remote address is a closed local port; no saved metrics
cp "$build/cravat-system.yml" "$pkg/conf/cravat-system.yml"
sed -n '/^base_modules:/,$p' "$pkg/cravat-system.template.yml" >> "$pkg/conf/cravat-system.yml"
# <site-packages>/cravat/metrics stays as an empty folder: cravat.constants creates it at import
rm -rf /opt/opencravat/metrics/* /opt/opencravat/workspace/* /opt/opencravat/jobs/* /opt/opencravat/logs/* "$pkg/metrics"/*
python -c "import cravat.constants"   # every folder it needs exists, so it writes nothing
chmod -R a+rX /opt/opencravat "$pkg"

# 3. checks: local addresses only, empty metrics folder, oc answers
python - "$pkg/conf/cravat-system.yml" <<'PY'
import sys, yaml
conf = yaml.safe_load(open(sys.argv[1]))
urls = {k: v for k, v in conf.items() if k.endswith("_url")}
assert urls and all(v == "http://127.0.0.1:9" for v in urls.values()), urls
assert conf["user_email_opt_out"] is True and conf["metrics_dir"] == "/opt/opencravat/metrics"
print("run-time configuration: every remote address is http://127.0.0.1:9")
PY
test -z "$(ls -A /opt/opencravat/metrics)" && test -z "$(ls -A "$pkg/metrics")"
oc version
ps --version
