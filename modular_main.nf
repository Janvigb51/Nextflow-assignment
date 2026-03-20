#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reads = '*_{1,2}.fq.gz'
params.outdir = 'outputs_modular/'
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
read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true).map { sample, reads -> tuple(sample, reads.collect {it.toAbsolutePath() }) }
adapter_ch = Channel.fromPath(params.adapters)
reference_ch = Channel.fromPath(params.reference, checkIfExists: true)
reference_prefix = file(params.reference).name


include { fastqc } from './modules/fastqc.nf'
include { trimmomatic } from './modules/trimmomatic.nf'
include { index_reference } from './modules/index.nf'
include { alignment } from './modules/alignment.nf'


// Run the workflow
workflow {
    read_pairs_ch.view()
    fastqc(read_pairs_ch)
    trimmomatic(read_pairs_ch, adapter_ch)
    indexed_ref_ch = index_reference(reference_ch)
    alignment(trimmomatic.out.trimmed_fq, indexed_ref_ch)
}
