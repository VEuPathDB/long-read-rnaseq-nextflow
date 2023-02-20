#!/usr/bin/env nextflow
nextflow.enable.dsl=2

sample_ch =   Channel
    .fromPath("$projectDir/data/Toxo/fastq/*fastq")
    //.splitFastq( by : 100, file:true  )

sample_ch.view()
paired_ch = Channel
    .fromFilePairs('/Users/saikouybah/Documents/RNA_Seq_Nexflow/data/Fastq/*_{1,2}.fq.gz', flat: true)
    //.splitFastq(by: 100_000, pe: true, file: true)
paired_ch.view()

params.reference1="/Users/saikouybah/Documents/RNA_Seq_Nexflow/data/Reference/pfal3D7_chr1.fasta"

process minimapMapping1{
    container = 'staphb/minimap2'
    
    input:
    path(reference)
    path(sample)

    output:
    //tuple val("${sample_base}"), path("*sam") 
    path("*sam") 

    script:
    sample_base = sample.getSimpleName()

    """
    minimap2 -ax map-ont --splice  -t 8 -2 --MD ${reference} ${sample} > ${sample}.sam
    """
}

process hisat2Index {
    container = 'veupathdb/shortreadaligner'
    input:
    path(reference)

    output:
    path 'genomeIndex*.ht2', emit: ht2_files
    val 'genomeIndex' , emit: genome_index_name

    script:

    """
    hisat2-build ${reference} genomeIndex
    """
}

process hisatMapping {
    container = 'veupathdb/shortreadaligner'
    
    input:
    path 'genomeIndex.*.ht2'
    val(genomeIndex)
    tuple val(sampleID), path(reads1), path(reads2)

    output:
    //tuple val("${sample_base}"), path("*sam") 
    path("*sam") 

    script:
    split_name = reads1.getBaseName()
    sample_base = reads1.getSimpleName()

    """
    hisat2 -x ${genomeIndex}  --max-intronlen 300 -1 ${reads1} -2 ${reads2} > ${split_name}.sam
    """

}

process sortSam {
    input:
    path(sam)

    output:
    tuple val("${sample_base}"), path("*sam")

    script:
    split_name = sam.getBaseName()
    sample_base = sam.getSimpleName()
    
    """
    samtools sort -n ${sam} -o ${split_name}_sorted.sam
    """
}
process mergesams {

    publishDir "${params.results}/sam", mode: 'copy'

    input:
    tuple val(sampleID), path("*.sam")

    output:
    path("${sampleID}.sam"), emit: sam
    val(sampleID), emit: sampleID

    script:

    """
    samtools merge -n -o ${sampleID}.sam *.sam
    """
}

workflow{
    sam =  minimapMapping1(params.reference, sample_ch)
    //sortedsam = sortSam(sam)
    //sortedsam.view()
    //samSet = sortedsam.groupTuple(sort: true)
    //samSet.view()
    //mergeSam = mergesams(samSet)
   //-------
    //index_ch = hisat2Index(params.reference1)
    //index_ch.ht2_files.view()
    //sample_ch.view()
    //paired_ch.view()
    //mapping = hisatMapping(index_ch.ht2_files, index_ch.genome_index_name, sample_ch)
    //sorted_sam = mapping.groupTuple()
    //sorted_sam.view()
    //mergesam = mergesams(sorted_sam)
    //mapping.name.view()
    //samSet = mapping.sam.collectFile()
    //mergesam = mergesams(mapping.sampleID, samSet)
    //mergesam.sam.view()
    //mergesam.sampleID.view()

}
