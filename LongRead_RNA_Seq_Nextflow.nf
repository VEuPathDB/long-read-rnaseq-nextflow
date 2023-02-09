#!/usr/bin/env nextflow
nextflow.enable.dsl=2


reads_ch = Channel.fromPath([params.reads + '/*.fastq', params.reads + '/*.fastq.gz', params.reads + '/*.fq.gz'])

process minimapMapping{
    
    input:
    path(reference)
    path(sample)

    output:
    path("*bam"), emit: bam
    val("${sample_base}"), emit: sampleID

    script:
    sample_base = sample.getSimpleName()

    """
    minimap2 -ax map-pb --MD "${reference}" "${sample}" > "${sample_base}".bam
    """
}


process transcriptClean {
    input:
    path(sam)
    path(reference)
    val(sample_base)

    output:
    path("${sample_base}_clean.sam")

    script:

    """
    python /usr/local/bin/TranscriptClean.py --sam "${sam}" --genome "${reference}" --outprefix "${sample_base}"
    """
}

process initiateDatabase {

    publishDir "${params.results}/database", mode: 'copy'

    input: 
    path(annotation)
    val(annot_name)
    val(build)
   

    output:
    path("talon.db")

    script:

    """
        talon_initialize_database --f "${annotation}" --a "${annot_name}" --g "${build}"  --o "talon"
    """
}

process talonLabelReads {

    input:
    path(sample)
    path(reference)
    val(sample_base)

    output:
    path("${sample_base}*.sam"), emit: samFiles
    val("${sample_base}"), emit: sample_base

    script:

    """
    talon_label_reads --f "${sample}" --g "${reference}" --t 1 --ar 20 --deleteTmp --o "${sample_base}"
    """

}

process genrateConfig {

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

    input:
    path(config)
    path(database)
    val(build)

    output:
    path("results_talon_read_annot.tsv"), emit: tsv_results
    path("results_QC.log")

    script:
    """   
    talon  --f "${config}" --db "${database}" --build "${build}" --o results
    """
}

process sampleList {
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

    input:
    path(database)
    path("results")

    output:
    path("results*")

    script:
    """
    talon_summarize --db "${database}" --v --o "results"
    """
}

process talonFilterTranscripts {

    input:
    path(database)
    path(datasets)
    val(annot_name) 

    output:
    path("filtered_transcripts.csv")

    script:
    """
    talon_filter_transcripts --db "${database}" --datasets "${datasets}" -a "${annot_name}" --maxFracA 0.5  --minCount 5  --minDatasets 2 --o filtered_transcripts.csv

    """
}

process transcriptAbundance{

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

    """
    talon_abundance --db "${database}" --whitelist "${wishList}"  -a "${annot_name}" --build "${build}" --o "results"
    """

}

process createGtf {

    publishDir "${params.results}/Gtf", mode: 'copy'

    input:
    path(database)
    val(annot_name)
    val(build)

    output:
    path("results*")

    """
    talon_create_GTF --db "${database}" -a "${annot_name}" --build "${build}" --o results
    """
}

process exctarctBysample{
    publishDir "${params.results}/counts", mode: 'copy'

    input:
    path(readCounts)

    output:
    path("*tsv")

    script:
    """
    subset_by_sample.py "${readCounts}"
    """
}
workflow{
    bam = minimapMapping(params.reference, reads_ch)
    cleanBam = transcriptClean(bam.bam,params.reference, bam.sampleID)

   database = initiateDatabase(params.referenceAnnotation, params.annotationName, params.build)
  
    label_reads = talonLabelReads(cleanBam, params.reference, bam.sampleID)
 
    samfiles = label_reads.samFiles.collect()
    samplesNames = label_reads.sample_base.collect()
    config = genrateConfig(samplesNames, params.build, params.platform, samfiles) 

    annotation = annotator(config.config_file, database, params.build)
 
    names_from_annotation = sampleList(annotation.tsv_results)
    talon_summary = talonSummarize(database, annotation.tsv_results)

    filtered = talonFilterTranscripts(database, names_from_annotation, params.annotationName) 

    abundance =  transcriptAbundance(database, filtered, params.annotationName, params.build, annotation.tsv_results)

    gtf = createGtf(database, params.annotationName, params.build)
    subsetCount = exctarctBysample(abundance)

}

