// Process alignment
process alignment {
    label "bwamem2"
    publishDir "${params.outdir}/aligned-reads-${sample}/", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(trimmed_reads)
    tuple path(reference_file), path(index_files)

    output:
    tuple val(sample), path("${sample}.sam")

    script:
    """
    bwa-mem2 mem ${reference_file} ${trimmed_reads[0]} ${trimmed_reads[1]} > ${sample}.sam
    """
}

