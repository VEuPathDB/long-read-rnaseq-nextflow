#!/usr/bin/env bash

set -euo pipefail

samtools merge -n -o ${sampleID}.sam *.sam

samtools view -S -b ${sampleID}.sam > tem.bam

samtools sort -n tem.bam -o ${sampleID}.bam
rm tem.bam

