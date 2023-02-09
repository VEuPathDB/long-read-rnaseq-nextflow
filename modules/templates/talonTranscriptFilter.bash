#!/usr/bin/env bash

set -euo pipefail

talon_filter_transcripts --db ${database} --datasets ${datasets} -a ${annot_name} --maxFracA 0.5  --minCount 5  --minDatasets 2 --o filtered_transcripts.csv