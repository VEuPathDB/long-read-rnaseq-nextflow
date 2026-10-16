#!/usr/bin/env nextflow

def fetchRunAccessions( tsv ) {
    def splitter = new nextflow.splitter.CsvSplitter().options( header:true, sep:'\t' )
    def reader = new BufferedReader( new FileReader( tsv ) )
    return splitter.target( reader ).list().collect { row -> row['run_accession'] }
}

//---------------------------------------
// include the RNA seq workflow
//---------------------------------------

include { longRna } from  './modules/longReadRnaSeq.nf'


//--------------------------------------
// Process the workflow
//-------------------------------------

workflow {
  if(!params.gtf) {
      error("Missing parameter params.gtf")
    }
  if(!params.fasta) {
      error("Missing parameter params.fasta")
    }

  if(!params.input) {
      error("Missing parameter params.input")
    }

  if(!params.platform) {
      error("Missing parameter params.platform")
    }

  if(!params.build) {
      error("Missing parameter params.build")
    }

  if(!params.annotationName) {
      error("Missing parameter params.annotationName")
    }
  if(!params.results) {
      error("Missing parameter params.results")
    }



  sampleRows = channel.fromPath(params.input + "/" + params.samplesheetFileName)
    .splitCsv( skip:1)

  sample_ch = sampleRows.map { row ->
    def fileName = file(row[1]);
    return [ [id: row[0] ], fileName ]
  }.map {
    return ([it[0], it[1].splitFastq(by : params.splitChunk, file:true )])
  }.transpose()

  longRna(sample_ch)
}
