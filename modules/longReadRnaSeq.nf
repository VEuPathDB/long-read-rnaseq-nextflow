#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process downloadSRA {
  container = 'veupathdb/bowtiemapping'

  input:
    val(sra)

  output:
    path("${sra}*")

  script:
    template 'fastqDump.bash'
}


process minimapMapping{
  container = 'staphb/minimap2'
    
  input:
    path(reference)
    path(sample)

  output:
    path("*sam") 

  script:
    sample_base = sample.getSimpleName()
    template 'minima2.bash'
}


process sortSam {
  container = 'veupathdb/shortreadaligner'

  input:
    path(sam)

  output:
    tuple val("${sample_base}"), path("*sam")

  script:
    split_name = sam.getBaseName()
    sample_base = sam.getSimpleName()
    template 'samSorting.bash'
}


process mergeSams {
  container = 'veupathdb/shortreadaligner'
    
  publishDir "${params.results}/bam", pattern: "*.bam",  mode: 'copy'

  input:
    tuple val(sampleID), path("*.sam")
    
  output:
    path("${sampleID}.sam"), emit: sam
    val(sampleID), emit: sampleID
    path("*bam"), emit: bam

  script:
    template 'samMerge.bash'  
}


process transcriptClean {
  container = 'sybah/longreadrnaseq'

  input:
    path(sam)
    path(reference)
    val(sample_base)

  output:
    path("${sample_base}_clean.sam")

  script:
    template 'transcriptClean.bash'
}


process initiateDatabase {
  container = 'sybah/longreadrnaseq'

  publishDir "$projectDir/data/database", mode: 'copy'

  input: 
    path(annotation)
    val(annot_name)
    val(build)
   
  output:
    path("*.db"), emit: db
    val(annot_name), emit: db_name

  script:
    template 'initDatabase.bash'
}


process talonLabelReads {
  container = 'sybah/longreadrnaseq'

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

process generateConfig {
    container = 'sybah/longreadrnaseq'

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


process annotator {
  container = 'sybah/longreadrnaseq'
    
  input:
    path(config)
    path(database)
    val(build)
    val(annot_name)

  output:
    path("results_talon_read_annot.tsv"), emit: tsv_results
    path("results_QC.log")

  script:
    template 'talonAnnotate.bash'
}


process sampleList {
  container = 'sybah/longreadrnaseq'

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
  container = 'sybah/longreadrnaseq'

  input:
    path(database)
    path("results")

  output:
    path("results*")

  script:
    template 'talonSummarise.bash'
}


process talonFilterTranscripts {
  container = 'sybah/longreadrnaseq'

  input:
    path(database)
    path(datasets)
    val(annot_name) 

  output:
    path("filtered_transcripts.csv")

  script:
    template 'talonTranscriptFilter.bash'
}


process transcriptAbundance{
  container = 'sybah/longreadrnaseq'

  publishDir "${params.results}/counts", mode: 'copy'

  input:
    path(database)
    path(wishList)
    val(annot_name)
    val(build)
    path("results*")

  output:
    path("results*")

  script:
    template 'talonAbundance.bash'
}


process transcriptAbundanceNoFilter{
  container = 'sybah/longreadrnaseq'

  publishDir "${params.results}/counts", mode: 'copy'

  input:
    path(database)
    val(annot_name)
    val(build)
    path("results*")

  output:
    path("results*")

  script:
    template 'talonAbundanceNoFilter.bash'
}


process createGtf {
  container = 'sybah/longreadrnaseq'

  publishDir "${params.results}/Gtf", mode: 'copy'

  input:
    path(annotOut)
    path(database)
    val(annot_name)
    val(build)

  output:
    path("*results*")

  script:
    template 'talonGtf.bash'    
}


process extractBysample{
  container = 'sybah/longreadrnaseq'
    
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

process indexGff {
    container = "veupathdb/proteintogenomealignment:latest"

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

    filtered = talonFilterTranscripts(params.database, namesFromAnnotation, params.annotationName) 

    abundanceNoFilter = transcriptAbundanceNoFilter(params.database, params.annotationName, params.build, annotation.tsv_results)
    abundanceFilter =  transcriptAbundance(params.database, filtered, params.annotationName, params.build, annotation.tsv_results)

    gtf = createGtf(annotation.tsv_results, params.database, params.annotationName, params.build)
    subsetCount = extractBysample(abundanceNoFilter, abundanceFilter )
    makeGff = convertGtfToGff(gtf)
    index = indexGff(makeGff)

}
