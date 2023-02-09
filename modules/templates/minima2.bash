#!/usr/bin/env bash

set -euo pipefail

minimap2 -ax map-pb --MD ${reference} ${sample} > ${sample_base}.bam