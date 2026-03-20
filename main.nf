#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reads = '*_{1,2}.fq.gz'
params.outdir = 'outputs/'
params.adapters = 'adapters.fa'
params.reference = 'LG12.fasta'

log.info """
      LIST OF PARAMETERS
================================
Reads            : ${params.reads}
Output-folder    : ${params.outdir}
Adapters         : ${params.adapters}
Reference-genome : ${params.reference}
"""

// Create channels
read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true).map { sample, reads -> tuple(sample, reads.collect { it.toAbsolutePath() }) }
adapter_ch = Channel.fromPath(params.adapters)
reference_ch = Channel.fromPath(params.reference, checkIfExists: true)
reference_prefix = file(params.reference).name

// Define fastqc process
process fastqc {
label "fastqc"
    publishDir "${params.outdir}/quality-control-${sample}/", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(reads)

    output:
    path("*_fastqc.{zip,html}")

    script:
    """
    fastqc ${reads}
    """
}

// Process trimmomatic
process trimmomatic {
label "trimmomatic"
    publishDir "${params.outdir}/trimmed-reads-${sample}/", mode: 'copy'

    input:
    tuple val(sample), path(reads)
    path adapters_file

    output:
    tuple val("${sample}"), path("${sample}*.trimmed.fq.gz"), emit: trimmed_fq
    tuple val("${sample}"), path("${sample}*.discarded.fq.gz"), emit: discarded_fq

    script:
    """
    trimmomatic PE -phred33 \
	${reads[0]} ${reads[1]} \
	${sample}_1.trimmed.fq.gz ${sample}_1.discarded.fq.gz \
	${sample}_2.trimmed.fq.gz ${sample}_2.discarded.fq.gz \
	ILLUMINACLIP:${adapters_file}:2:30:10
    """
}

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


// Run the workflow
workflow {
    read_pairs_ch.view()
    fastqc(read_pairs_ch)
    trimmomatic(read_pairs_ch, adapter_ch)
    indexed_ref_ch = index_reference(reference_ch)
    alignment(trimmomatic.out.trimmed_fq, indexed_ref_ch)
}

