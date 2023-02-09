#!/usr/bin/env bash

set -euo pipefail

python /usr/local/bin/TranscriptClean.py --sam ${sam} --genome ${reference} --outprefix ${sample_base}