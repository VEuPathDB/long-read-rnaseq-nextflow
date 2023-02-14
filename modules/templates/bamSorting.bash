#!/usr/bin/env bash

set -euo pipefail

samtools sort ${bam} -o ${sample_base}.bam