#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
* This process download the Fastq sequence files from the sequence read archive

*/
process downloadSRA {
  container = 'veupathdb/bowtiemapping:v1.0.0'

  input:
    val(sra)

  output:
    path("${sra}*")

  script:
    template 'fastqDump.bash'
}

/*
This below process map the long read RNA-Seq data to the reference genome using minimap2 

@reference is the refence genome in fasta format
@sample in the sample to be mapped

For efficiency in mapping, each fastq file is split in to smaller chuck, mapp, coordinate sort and merge.

The out put is a sam file
*/

process minimapMapping{
  container = 'staphb/minimap2:2.28'
    
  input:
    path(reference)
    path(sample)

  output:
    path("*sam") 

  script:
    sample_base = sample.getSimpleName()
    template 'minima2.bash'
}
/*
This process sort the alignment file (sam) by cooridinates

@sam is the sam file generated from the mapping step above

Output a coordinate sorted sam file.
*/

process sortSam {
  container = 'veupathdb/shortreadaligner:v1.0.0'

  input:
    path(sam)

  output:
    tuple val("${sample_base}"), path("*sam")

  script:
    split_name = sam.getBaseName()
    sample_base = sam.getSimpleName()
    template 'samSorting.bash'
}

/*

This process merge the coordinate sorted sam files
@sampleID tuple of sample IDs of sam file from the same split above

Output is bam file of merge sam files. 

*/

process mergeSams {
  container = 'veupathdb/shortreadaligner:v1.0.0'
    
  publishDir "${params.results}/bam", pattern: "*.bam*",  mode: 'copy'

  input:
    tuple val(sampleID), path("*.sam")
    
  output:
    path("${sampleID}.sam"), emit: sam
    val(sampleID), emit: sampleID
    path("*bam"), emit: bam
    path("*bam.bai"), emit: bam_bai

  script:
    template 'samMerge.bash'  
}

/*
This process run TranscriptClean to fix non-canonical jubctions
@sam, coordinate sorted sam file
@reference, reference genome in fasta format
@sample_base, sample base name

Output is sam file with corrected non-canonical junctions

*/

process transcriptClean {
  container = 'veupathdb/longreadrnaseq:v1.0.0'

  input:
    path(sam)
    path(reference)
    val(sample_base)

  output:
    path("${sample_base}_clean.sam")

  script:
    template 'transcriptClean.bash'
}

/*
This process initialise the TALON database using the current available annotation. 

*/
process initiateDatabase {
  container = 'veupathdb/longreadrnaseq:v1.0.0'

  publishDir "${params.databaseDir}", mode: 'copy'

  input: 
    path(annotation)
    val(annotationName)
    val(build)
   
  output:
    path("*.db"), emit: db
    val(annotationName), emit: db_name

  script:
    template 'initDatabase.bash'
}

/*
This process label the reads with potential priming sites
*/
process talonLabelReads {
  container = 'veupathdb/longreadrnaseq:v1.0.0'

  input:
    path(sample)
    path(reference)
    val(sample_base)

  output:
    path("${sample_base}*.sam"), emit: samFiles
    val("${sample_base}"), emit: sample_base

  script:
    template 'readLabel.bash'
}

/*
This process generate the TALON configuration file
*/

process generateConfig {
    container = 'veupathdb/longreadrnaseq:v1.0.0'

  input:
    val(samID)
    val(build)
    val(seqPlatform)
    val(sampleName)

  output:
    path("config.txt"), emit: config_file
    path("Sample_names.txt"), emit: sampleNames

  script:
    """
    config.py "${samID}" "${build}"  "${seqPlatform}"  "${sampleName}"
    """
}

/*

This process annotate the transcripts 
*/

process annotator {
  container = 'veupathdb/longreadrnaseq:v1.0.0'
    
  input:
    path(config)
    path(database)
    val(build)
    val(annotationName)

  output:
    path("results_talon_read_annot.tsv"), emit: tsv_results
    path("results_QC.log")

  script:
    template 'talonAnnotate.bash'
}

/*
This process generate the sample sample list from the annotation database

*/
process sampleList {
  container = 'veupathdb/longreadrnaseq:v1.0.0'

  input:
    path(annotation)

  output:
    path("Sample_names.txt")

  script:
    """
    gene_names.py "${annotation}"
    """
}


process talonSummarize {
  container = 'veupathdb/longreadrnaseq:v1.0.0'

  input:
    path(database)
    path("results")

  output:
    path("results*")

  script:
    template 'talonSummarise.bash'
}

/*
Apply filter to TALON transcript using these talon default setting maxFracA = 0.5, minCount = 5, minDatasets = 2

*/

process talonFilterTranscripts {
  container = 'veupathdb/longreadrnaseq:v1.0.0'

  input:
    path(database)
    path(datasets)
    val(annotationName)
    val(maxFracA)
    val(minCount)
    val(minDatasets)

  output:
    path("filtered_transcripts.csv")

  script:
    template 'talonTranscriptFilter.bash'
}

/*
Determine transcript abudance for individual transcripts using TALON default filter from the above process and put them a matrix
*/
process transcriptAbundance{
  container = 'veupathdb/longreadrnaseq:v1.0.0'

  publishDir "${params.results}/counts", mode: 'copy'

  input:
    path(database)
    path(wishList)
    val(annotationName)
    val(build)
    path("results*")

  output:
    path("results*")

  script:
    template 'talonAbundance.bash'
}

/*
Filter transcript abudance for individual transcripts without a filter and put them a matrix
*/

process transcriptAbundanceNoFilter{
  container = 'veupathdb/longreadrnaseq:v1.0.0'

  publishDir "${params.results}/counts", mode: 'copy'

  input:
    path(database)
    val(annotationName)
    val(build)
    path("results*")

  output:
    path("results*")

  script:
    template 'talonAbundanceNoFilter.bash'
}

/*
Generate an annotation file (Gtf) based on the gene model identified by talon
*/


process createGtf {
  container = 'veupathdb/longreadrnaseq:v1.0.0'

  publishDir "${params.results}/Gtf", mode: 'copy'

  input:
    path(annotOut)
    path(database)
    val(annotationName)
    val(build)

  output:
    path("*results*")

  script:
    template 'talonGtf.bash'    
}

/*
Extract results from individual samples from the expression matrix genetated by TALON
*/
process extractBysample{
  container = 'veupathdb/longreadrnaseq:v1.0.0'
    
  publishDir "${params.results}/counts", mode: 'copy'

  input:
    path(unFilteredCounts)
    path(filteredCounts)

  output:
    path("*tsv")

  script:
    """
    subset_by_sample.py "${unFilteredCounts}" "${filteredCounts}"
    """
}

/*
Convert the TALON generated Gtf into Gff
*/

process convertGtfToGff {
  container = 'quay.io/biocontainers/agat:0.8.0--pl5262hdfd78af_0'

  publishDir "${params.results}/Gtf", mode: 'copy'

  input:
    path(gtf)

  output:
    path("*gff")

    script:
    template 'gtfToGff.bash'
    
}

/*
Process indix the final gtf file

*/
process indexGff {
    container = "veupathdb/proteintogenomealignment:v1.0.0"

    publishDir "${params.results}/Gtf", mode: 'copy'

    input:
    path(gff)

    output:
    path("*gff*")
    
    script:
    template 'indexGff.bash'
}


workflow longRna {
  take: 
    sample_ch

  main:
    if (params.local) {
      sam =  minimapMapping(params.reference, sample_ch)
    } else {
      sample = downloadSRA(sample_ch)
      .splitFastq( by : params.splitChunk, file:true )  
      sam =  minimapMapping(params.reference, sample)
    }
       
    sortedsam = sortSam(sam)
    samSet = sortedsam.groupTuple(sort: true)

    mergeSam = mergeSams(samSet)
    cleanSam = transcriptClean(mergeSam.sam,params.reference, mergeSam.sampleID)

    database = initiateDatabase(params.referenceAnnotation, params.annotationName, params.build)

    labelReads = talonLabelReads(cleanSam, params.reference, mergeSam.sampleID)
    
    samfiles = labelReads.samFiles.collect()
    samplesNames = labelReads.sample_base.collect()
    config = generateConfig(samplesNames, params.build, params.platform, samfiles) 

    annotation = annotator(config.config_file, params.database, params.build, database.db_name)

    namesFromAnnotation = sampleList(annotation.tsv_results)
    talonSummary = talonSummarize(params.database, annotation.tsv_results)

    filtered = talonFilterTranscripts(params.database, namesFromAnnotation, params.annotationName, params.maxFracA, params.minCount, params.minDatasets) 

    abundanceNoFilter = transcriptAbundanceNoFilter(params.database, params.annotationName, params.build, annotation.tsv_results)
    abundanceFilter =  transcriptAbundance(params.database, filtered, params.annotationName, params.build, annotation.tsv_results)

    gtf = createGtf(annotation.tsv_results, params.database, params.annotationName, params.build)
    subsetCount = extractBysample(abundanceNoFilter, abundanceFilter )
    makeGff = convertGtfToGff(gtf)
    index = indexGff(makeGff)

}
