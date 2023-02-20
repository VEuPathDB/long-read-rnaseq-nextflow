#!/usr/bin/env bash

set -euo pipefail

samtools merge -n -o ${sampleID}.sam *.sam
