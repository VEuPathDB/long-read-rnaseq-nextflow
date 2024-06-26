# <p align=center>Long Read RNA-Seq transcript annotation Workflow</p>
This workflow identify and quantify know and novel genes/transcript isoforms using long read RNA-seq data. The [TALON](https://github.com/mortazavilab/TALON) pipeline is used for the transcript identification and quantification. The long read are mapped to a reference genome using minimap2 and the sam file sorted using Samtools. After sorting [TranscriptClean](https://github.com/mortazavilab/TranscriptClean) is used to fix noncanonical junctions. 
TALON is then used the label the read to identify potential internal priming. A database of existing gene models are then generated based on the annotation to identify existing gene model and potential novel models. After identification of gene models, TALON quantify the expression level of each gene transcript. More details on the talon method can be found [here](https://www.biorxiv.org/content/10.1101/672931v2.full)


**<p align=left>Get Started</p>**
To run the work the following dependencies need to be install
* Docker
> `https://docs.docker.com/engine/install/`
* Nextflow
> `curl https://get.nextflow.io | bash`

* The pull the git hub repo using the following command
> `git pull https://github.com/VEuPathDB/bulk-rnaseq-nextflow.git`

* Alternatively the workflow can be run directly using nextflow which pull down the repo. 
> `nextflow run VEuPathDB/bulk-rnaseq-nextflow -with-trace -c  <config_file> -r main`


<br />


**<p align=left>Input Data</p>**
Example of input data can be found in the `data` directory. The following files are require
* Long read Fastq files or a csv file containing the list SRA accession numbers to be analyzed.
* A reference for the organism being analyzed
* A nextflow config file `nextflow.config`.

**<p align=left>Input Data</p>**
The following output files are generated and can be found in the result folder specify in the config file
* Un-filtered transcripts counts
* Filtered transcripts counts
* GTF annotation file generated from un-filters annotations
* GTF file generated from filtered annotation. 
* A bam file for each sample analyzed


***<p align=center>Nextflow workflow diagram</p>***  
```mermaid
flowchart TD
    p0((Channel.fromPath))
    p1([splitFastq])
    p2(( ))
    p3[longRna:minimapMapping]
    p4[longRna:sortSam]
    p5([groupTuple])
    p6[longRna:mergeSams]
    p8(( ))
    p9[longRna:transcriptClean]
    p10(( ))
    p11(( ))
    p12(( ))
    p13[longRna:initiateDatabase]
    p15(( ))
    p16[longRna:talonLabelReads]
    p17([collect])
    p18([collect])
    p19(( ))
    p20(( ))
    p21[longRna:generateConfig]
    p23(( ))
    p24(( ))
    p25[longRna:annotator]
    p27[longRna:sampleList]
    p28(( ))
    p29[longRna:talonSummarize]
    p30(( ))
    p31(( ))
    p32(( ))
    p33[longRna:talonFilterTranscripts]
    p34(( ))
    p35(( ))
    p36(( ))
    p37[longRna:transcriptAbundanceNoFilter]
    p38(( ))
    p39(( ))
    p40(( ))
    p41[longRna:transcriptAbundance]
    p42(( ))
    p43(( ))
    p44(( ))
    p45[longRna:createGtf]
    p46[longRna:extractBysample]
    p47(( ))
    p48[longRna:convertGtfToGff]
    p49[longRna:indexGff]
    p50(( ))
    p0 --> p1
    p1 -->|sample_ch| p3
    p2 -->|reference| p3
    p3 --> p4
    p4 --> p5
    p5 -->|samSet| p6
    p6 --> p9
    p6 -->|sampleID| p9
    p8 -->|reference| p9
    p9 --> p16
    p10 -->|annotation| p13
    p11 -->|annot_name| p13
    p12 -->|build| p13
    p13 -->|annot_name| p25
    p15 -->|reference| p16
    p6 -->|sample_base| p16
    p16 --> p17
    p16 --> p18
    p17 -->|samfiles| p21
    p18 -->|samplesNames| p21
    p19 -->|build| p21
    p20 -->|seqPlatform| p21
    p21 --> p25
    p23 -->|database| p25
    p24 -->|build| p25
    p25 --> p27
    p27 --> p33
    p28 -->|database| p29
    p25 -->|results| p29
    p29 --> p30
    p31 -->|database| p33
    p32 -->|annot_name| p33
    p33 --> p41
    p34 -->|database| p37
    p35 -->|annot_name| p37
    p36 -->|build| p37
    p25 -->|results*| p37
    p37 --> p46
    p38 -->|database| p41
    p39 -->|annot_name| p41
    p40 -->|build| p41
    p25 -->|results*| p41
    p41 --> p46
    p25 -->|annotOut| p45
    p42 -->|database| p45
    p43 -->|annot_name| p45
    p44 -->|build| p45
    p45 --> p48
    p46 --> p47
    p48 --> p49
    p49 --> p50
```