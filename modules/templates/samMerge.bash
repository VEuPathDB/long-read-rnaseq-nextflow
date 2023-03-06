#!/usr/bin/env bash

set -euo pipefail

samtools merge -f  ${sampleID}_temp.sam *.sam

samtools sort ${sampleID}_temp.sam -o ${sampleID}.sam

samtools view -bS ${sampleID}.sam > ${sampleID}.bam

rm ${sampleID}_temp.sam