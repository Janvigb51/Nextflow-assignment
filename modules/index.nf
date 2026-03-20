//Process index
process index_reference {
    label "bwamem2"
    publishDir "${params.outdir}/reference-index/", mode: 'copy', overwrite: true

    input:
    path reference_file

    output:
    tuple path(reference_file), path("${reference_file}.*")

    script:
    """
    bwa-mem2 index ${reference_file}
    """
}
