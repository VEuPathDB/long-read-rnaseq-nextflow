#!/usr/bin/env nextflow
nextflow.enable.dsl=2


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
  tuple val(meta), path(sample)

  output:
  tuple val(meta), path("minimap.sam")

  script:
  """
  minimap2 -ax splice -k14 -uf -2 -G 5000  ${reference} ${sample} > minimap.sam
  """

}
/*
This process sort the alignment file (sam) by cooridinates

@sam is the sam file generated from the mapping step above

Output a coordinate sorted sam file.
*/

process sortSam {
  container = 'quay.io/biocontainers/samtools:1.20--h50ea8bc_0'

  input:
  tuple val(meta), path(sam)

  output:
  tuple val(meta), path("sorted.sam")

  script:
  """
  samtools sort ${sam} -o sorted.sam
  """
}

/*

This process merge the coordinate sorted sam files
@sampleID tuple of sample IDs of sam file from the same split above

Output is bam file of merge sam files. 

*/

process mergeSams {
  container = 'quay.io/biocontainers/samtools:1.20--h50ea8bc_0'

  input:
  tuple val(meta), path("*.sam")  
    
  output:
  tuple val(meta), path("merged_sorted.sam")
  
  script:
  """
  samtools merge -f merged.sam *.sam
  samtools sort merged.sam -o merged_sorted.sam
  """
}


process bam {
  container = 'quay.io/biocontainers/samtools:1.20--h50ea8bc_0'

  publishDir "${params.results}/bam", pattern: "*.bam*",  mode: 'copy'
  
  input:
  tuple val(meta), path(sam)

  output:
  path("${meta.id}.bam")
  path("${meta.id}.bam.bai")

  script:
  """
  samtools view -bS $sam > ${meta.id}.bam
  samtools index ${meta.id}.bam
  """
}


/*
This process run TranscriptClean to fix non-canonical jubctions
@sam, coordinate sorted sam file
@reference, reference genome in fasta format
@sample_base, sample base name

Output is sam file with corrected non-canonical junctions

*/

process transcriptClean {
  container = 'veupathdb/longreadrnaseq:1.0.0'

  input:
  tuple val(meta), path(sam)
  path(reference)
  
  output:
  tuple val(meta), path("${meta.id}_clean.sam")

  script:
  """
  python /usr/local/bin/TranscriptClean.py --sam ${sam} \
    --genome ${reference} \
    --outprefix ${meta.id}
  """

}

/*
This process initialise the TALON database using the current available annotation. 

*/
process initiateDatabase {
  container = 'veupathdb/longreadrnaseq:1.0.0'

  input:
  path(annotation)
  val(annotationName)
  val(build)
   
  output:
  path("${build}.db")

  script:
  """
  talon_initialize_database --f ${annotation} \
    --a ${annotationName} \
    --g ${build}  \
    --o ${build}
  """
}

/*
This process label the reads with potential priming sites
*/
process talonLabelReads {
  container = 'veupathdb/longreadrnaseq:1.0.0'

  input:
  tuple val(meta), path(sam)
  path(reference)

  output:
  tuple val(meta), path("${meta.id}_labeled.sam")
  
  script:
  """
  talon_label_reads --f ${sam} \
                  --g ${reference} \
                  --t 1 \
                  --ar 20 \
                  --deleteTmp \
                  --o ${meta.id}
  """
}

/*
This process generate the TALON configuration file
*/

process generateConfig {
    container = 'veupathdb/longreadrnaseq:1.0.0'

  input:
    tuple val(meta), path(sam)
    val(build)
    val(platform)

  output:
    path("config.txt")

  script:
  """
  samFullPath=\$(readlink -f $sam)
  echo ${meta.id},$build,$platform,\$samFullPath >config.txt
  """
}

/*

This process annotate the transcripts 
*/

process annotator {
  container = 'veupathdb/longreadrnaseq:1.0.0'
    
  input:
  path(configFile)
  path(database)
  val(build)

  output:
  path("results_talon_read_annot.tsv"), emit: results
  path("${database}_modified.db"), emit: database

  script:
  """
  cp $database ${database}_modified.db
  talon --f ${configFile} \
    --db ${database}_modified.db \
    --build ${build} \
    --o results
  """

}

/*
This process generate the sample sample list from the annotation database

*/
process sampleList {
  container = 'veupathdb/longreadrnaseq:1.0.0'

  input:
    path(annotation)

  output:
    path("Sample_names.txt")

  script:
    """
    gene_names.py "${annotation}"
    """
}


/*
Apply filter to TALON transcript using these talon default setting maxFracA = 0.5, minCount = 5, minDatasets = 2

*/

process talonFilterTranscripts {
  container = 'veupathdb/longreadrnaseq:1.0.0'

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
  """
  talon_filter_transcripts --db ${database}\
    --datasets ${datasets} \
    -a ${annotationName} \
    --maxFracA ${maxFracA}  \
    --minCount ${minCount}  \
    --minDatasets ${minDatasets} \
    --o filtered_transcripts.csv
  """
}

/*
Determine transcript abudance for individual transcripts using TALON default filter from the above process and put them a matrix
*/
process transcriptAbundance{
  container = 'veupathdb/longreadrnaseq:1.0.0'

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
  """
  talon_abundance --db ${database} \
    --whitelist ${wishList}  \
    -a ${annotationName} \
    --build ${build} \
    --o "results"
  """

}

/*
Filter transcript abudance for individual transcripts without a filter and put them a matrix
*/

process transcriptAbundanceNoFilter{
  container = 'veupathdb/longreadrnaseq:1.0.0'

  publishDir "${params.results}/counts", mode: 'copy'

  input:
  path(database)
  val(annotationName)
  val(build)
  path("results*")

  output:
  path("results*")

  script:
  """
  talon_abundance --db ${database}  \
    -a ${annotationName} \
    --build ${build} \
    --o "results_no_filter"
  """

}

/*
Generate an annotation file (Gtf) based on the gene model identified by talon
*/
process createGtf {
  container = 'veupathdb/longreadrnaseq:1.0.0'

  publishDir "${params.results}/Gtf", mode: 'copy'

  input:
    path(database)
    val(annotationName)
    val(build)

  output:
    path("*results*")

  script:
  """
  talon_create_GTF --db ${database} \
    --observed \
    -a ${annotationName} \
    --build ${build} \
    --o ${build}_results
  """
}

/*
Extract results from individual samples from the expression matrix genetated by TALON
*/
process extractBysample{
  container = 'veupathdb/longreadrnaseq:1.0.0'
    
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
  """
  agat_convert_sp_gxf2gxf.pl -g ${gtf} \
    -o ${params.build}.gff
  """

}

/*
Process indix the final gtf file

*/
process indexGff {
  container = 'biocontainers/tabix:v1.9-11-deb_cv1'

  publishDir "${params.results}/Gtf", mode: 'copy'

  input:
  path(gff)

  output:
  path("*gff*")
    
  script:
  """
  sort -k1,1 -k4,4n ${gff} > ${params.build}_sorted.gff
  bgzip ${params.build}_sorted.gff
  tabix -p gff ${params.build}_sorted.gff.gz
  """
}


/*
Rewrite the fasta file with a clean defline
*/
process createGtf {
  container = 'veupathdb/longreadrnaseq:1.0.0'

  input:
    path(fasta)

  output:
    path("cleanedGenomic.fasta")

  script:
  """
  perl -e 'while(<>){chomp; if(/(^>\S+)/){print \$1 . "\n"}else {print \$_ . "\n";}}' $fasta >cleanedGenomic.fasta
  """
}



workflow longRna {
  take: 
    sample_ch

  main:


  cleanFasta = cleanFastaDefline(params.fasta)
  
  sam =  minimapMapping(cleanFasta, sample_ch)

  sortedsam = sortSam(sam)
  samSet = sortedsam.groupTuple()

  mergeSam = mergeSams(samSet)
  cleanSam = transcriptClean(mergeSam, cleanFasta)

  bam(cleanSam)
  
  initDatabase = initiateDatabase(params.gtf, params.annotationName, params.build)

  labelReads = talonLabelReads(cleanSam, cleanFasta)

  config = generateConfig(labelReads, params.build, params.platform)

  annotation = annotator(config.collectFile(), initDatabase, params.build)

  namesFromAnnotation = sampleList(annotation.results)

  filtered = talonFilterTranscripts(annotation.database, namesFromAnnotation, params.annotationName, params.maxFracA, params.minCount, params.minDatasets)

  abundanceNoFilter = transcriptAbundanceNoFilter(annotation.database, params.annotationName, params.build, annotation.results)
  abundanceFilter =  transcriptAbundance(annotation.database, filtered, params.annotationName, params.build, annotation.results)

  gtf = createGtf(annotation.database, params.annotationName, params.build)
  subsetCount = extractBysample(abundanceNoFilter, abundanceFilter )
  makeGff = convertGtfToGff(gtf)
  index = indexGff(makeGff)

}


