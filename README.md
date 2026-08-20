# long-read-rnaseq-nextflow

Nextflow pipeline that identifies and quantifies known and novel gene/transcript isoforms from long-read RNA-seq data for VEuPathDB genome resources.

## Overview

Long reads are split into chunks and mapped to a reference genome with `minimap2`, coordinate-sorted and merged per sample with `samtools`, then corrected for noncanonical splice junctions with `TranscriptClean`. [TALON](https://github.com/mortazavilab/TALON) labels reads for potential internal priming, initializes a transcript database from the existing reference annotation, and annotates each sample's reads against it to identify known and novel gene models. TALON's abundance tools then quantify each transcript both with and without TALON's default noise filters (`maxFracA`, `minCount`, `minDatasets`), and the resulting gene models are written out as a GTF, converted to a sorted/indexed GFF3, and split back out per sample. Background on the TALON method is available [here](https://www.biorxiv.org/content/10.1101/672931v2.full).

This pipeline produces the transcript models and per-sample/aggregate expression matrices used to annotate and display long-read RNA-seq datasets on VEuPathDB genome sites.

## Requirements

- Nextflow
- Docker (default) or Singularity/Apptainer — see `conf/docker.config`, `conf/singularity.config`
- Processes run in `staphb/minimap2`, `quay.io/biocontainers/samtools`, `veupathdb/longreadrnaseq` (TALON, TranscriptClean, and pipeline helper scripts), `quay.io/biocontainers/agat` (GTF-to-GFF3 conversion), and `biocontainers/tabix` (GFF indexing)

## Usage

```
nextflow run VEuPathDB/long-read-rnaseq-nextflow -r main \
  --input /path/to/fastq_dir \
  --samplesheetFileName samplesheet.csv \
  --fasta /path/to/reference.fa \
  --gtf /path/to/reference.gtf \
  --build MyOrganismDB \
  --annotationName MyOrganismDB \
  --platform Nanopore \
  --results /path/to/results \
  -resume -C site.config
```

`samplesheet.csv` (found via `input`/`samplesheetFileName`) is a CSV with a header row followed by one line per sample: sample ID in the first column, path to that sample's FASTQ file in the second.

The pipeline has a single entry point (the default, unnamed `workflow`); there are no named `-entry` targets.

## Key parameters

- `input` — directory containing the input FASTQ files and the sample sheet
- `samplesheetFileName` — name of the sample sheet CSV within `input` (sample ID, FASTQ path)
- `fasta` — reference genome FASTA (cleaned of extra defline text before mapping)
- `gtf` — reference annotation GTF, used to initialize the TALON database
- `build` — genome build/assembly name recorded in the TALON database and used to name output files
- `annotationName` — annotation name recorded in the TALON database
- `platform` — sequencing platform label recorded in the TALON config (e.g. `Nanopore`)
- `splitChunk` — number of reads per FASTQ chunk when splitting each sample for parallel mapping
- `maxFracA` / `minCount` / `minDatasets` — TALON transcript filter thresholds (fraction of A's at the 3' end, minimum read count, minimum number of datasets a transcript must appear in) applied when producing the filtered abundance/whitelist
- `results` — output directory for BAMs, GTF/GFF, and count files

## Output

Written under `results`:

- `bam/` — coordinate-sorted, indexed BAM per sample
- `Gtf/` — TALON-generated GTF of observed gene models, its AGAT-converted GFF3, and a sorted, `bgzip`-compressed, `tabix`-indexed version of that GFF3
- `counts/` — unfiltered and TALON-filtered transcript abundance matrices, plus per-sample abundance files split from the aggregate matrices

Filtered outputs apply TALON's `maxFracA`/`minCount`/`minDatasets` thresholds to distinguish high-confidence transcripts from potential internal-priming or low-support artifacts; unfiltered outputs retain the full set of observed transcripts.
