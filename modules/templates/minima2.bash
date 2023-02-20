#!/usr/bin/env bash

set -euo pipefail

minimap2 -ax map-ont --splice  -t 8 -2 --MD ${reference} ${sample} > ${sample}.sam
