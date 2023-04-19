#!/usr/bin/env nextflow
import nextflow.splitter.CsvSplitter

def fetchRunAccessions( tsv ) {
    def splitter = new CsvSplitter().options( header:true, sep:'\t' )
    def reader = new BufferedReader( new FileReader( tsv ) )
    splitter.parseHeader( reader )
    List<String> run_accessions = []
    Map<String,String> row
    while( row = splitter.fetchRecord( reader ) ) {
       run_accessions.add( row['run_accession'] )
    }
    return run_accessions
}

//---------------------------------------
// include the RNA seq workflow
//---------------------------------------

include { longRna } from  './modules/longReadRnaSeq.nf'


//---------------------------------------------------------------
// Param Checking 
//---------------------------------------------------------------


if(!params.referenceAnnotation) {
    throw new Exception("Missing parameter params.referenceAnnotation")
  }
if(!params.reference) {
    throw new Exception("Missing parameter params.reference")
  }

if(!params.reads && !params.sraAccession) {
    throw new Exception("Missing parameter params.reads and parameter params.sraAccession")
  }

if(!params.platform) {
    throw new Exception("Missing parameter params.platform")
  }

if(!params.build) {
    throw new Exception("Missing parameter params.build")
  }

if(!params.annotationName) {
    throw new Exception("Missing parameter params.annotationName")
  }
if(!params.results) {
    throw new Exception("Missing parameter params.results")
  }

if (params.local){
    sample_ch =   Channel
      .fromPath([params.reads + '/*.fastq', params.reads + '/*.fastq.gz', params.reads + '/*.fq.gz'])
      .splitFastq( by : params.splitChunk, file:true  )
} else {
    input = fetchRunAccessions(params.sraAccession)
    sample_ch = Channel.fromList(input)
}
//--------------------------------------
// Process the workflow
//-------------------------------------

workflow {
    longRna(sample_ch)
}