#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process minimapMapping{
    
    input:
    path(reference)
    path(sample)

    output:
    path("*bam"), emit: bam
    val("${sample_base}"), emit: sampleID

    script:
    sample_base = sample.getSimpleName()

    template 'minima2.bash'
}


process transcriptClean {
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

    publishDir "${params.results}/database", mode: 'copy'

    input: 
    path(annotation)
    val(annot_name)
    val(build)
   

    output:
    path("talon.db")

    script:
    template 'initDatabase.bash'
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
    template 'readLabel.bash'

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
    template 'talonAnnotate.bash'
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
    template 'talonSummarise.bash'
}

process talonFilterTranscripts {

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

    publishDir "${params.results}/Gtf", mode: 'copy'

    input:
    path(database)
    val(annot_name)
    val(build)

    output:
    path("results*")

    script:
    template 'talonGtf.bash'
    
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
workflow longRna {
    take: 
        reads_ch

    main:
        bam = minimapMapping(params.reference, reads_ch)
        cleanBam = transcriptClean(bam.bam,params.reference, bam.sampleID)

        database = initiateDatabase(params.referenceAnnotation, params.annotationName, params.build)
    
        labelReads = talonLabelReads(cleanBam, params.reference, bam.sampleID)
    
        samfiles = labelReads.samFiles.collect()
        samplesNames = labelReads.sample_base.collect()
        config = genrateConfig(samplesNames, params.build, params.platform, samfiles) 

        annotation = annotator(config.config_file, database, params.build)
        annotation.tsv_results.view()

        namesFromAnnotation = sampleList(annotation.tsv_results)
        talonSummary = talonSummarize(database, annotation.tsv_results)

        filtered = talonFilterTranscripts(database, namesFromAnnotation, params.annotationName) 

        abundance =  transcriptAbundance(database, filtered, params.annotationName, params.build, annotation.tsv_results)
        abundanceNoFilter = transcriptAbundanceNoFilter(database, params.annotationName, params.build, annotation.tsv_results)

        gtf = createGtf(database, params.annotationName, params.build)
        subsetCount = exctarctBysample(abundance)

}

