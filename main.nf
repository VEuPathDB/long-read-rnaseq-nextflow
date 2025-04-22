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


if(!params.gtf) {
    throw new Exception("Missing parameter params.gtf")
  }
if(!params.fasta) {
    throw new Exception("Missing parameter params.fasta")
  }

if(!params.input) {
    throw new Exception("Missing parameter params.input")
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



sampleRows = Channel.fromPath(params.input + "/" + params.samplesheetFileName)
  .splitCsv( skip:1)

sample_ch = sampleRows.map { row ->
  fileName = file(params.input + "/" + row[1]);
  return [ [id: row[0] ], fileName ]
}.map {
  return ([it[0], it[1].splitFastq(by : params.splitChunk, file:true )])
}.transpose()

//--------------------------------------
// Process the workflow
//-------------------------------------

workflow {
  longRna(sample_ch)
}
