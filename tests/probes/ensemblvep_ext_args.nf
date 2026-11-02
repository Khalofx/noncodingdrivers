//
// Test probe for tests/vep_ext_args.nf.test, never run by the pipeline: a native process with the
// name of the nf-core module, so the ENSEMBLVEP_VEP selector of conf/modules.config applies to it.
// It writes the ext.args string that the module receives for the same inputs, evaluated exactly
// as for the real task (the closure reads params and the extra_files input).
//
process ENSEMBLVEP_VEP {
    tag "${meta.id}"

    input:
    tuple val(meta), path(vcf), path(custom_extra_files)
    path extra_files

    output:
    tuple val(meta), path('ext_args.txt'), emit: args

    exec:
    task.workDir.resolve('ext_args.txt').text = task.ext.args
}
