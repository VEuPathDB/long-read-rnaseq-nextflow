#!/usr/bin/env bash

set -euo pipefail

talon_filter_transcripts --db ${database}\
                         --datasets ${datasets} \
                         -a ${annotationName} \
                         --maxFracA ${maxFracA}  \
                         --minCount ${minCount}  \
                         --minDatasets ${minDatasets} \
                         --o filtered_transcripts.csv