#!/usr/bin/env bash

set -euo pipefail

minimap2 -ax splice -k14 -uf -G 5000  ${reference} ${sample} > ${sample}.sam
