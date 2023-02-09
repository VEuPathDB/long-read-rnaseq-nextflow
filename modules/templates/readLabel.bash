#!/usr/bin/env bash

set -euo pipefail

talon_label_reads --f ${sample} --g ${reference} --t 1 --ar 20 --deleteTmp --o ${sample_base}