#!/usr/bin/env bash

set -euo pipefail

minimap2 -ax map-pb -t 8 -2 --MD ${reference} ${sample} > ${sample_base}.bam