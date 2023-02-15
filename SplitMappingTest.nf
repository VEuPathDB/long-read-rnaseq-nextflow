#!/usr/bin/env nextflow
nextflow.enable.dsl=2
sample_ch =   Channel
    .fromPath("$projectDir/data/Toxo/fastq/*fastq")
    .splitFastq( record: true, by : 10  )
    .view { record -> record.readHeader }