#!/usr/bin/env bash

set -euo pipefail

talon_abundance --db ${database}  -a ${annotationName} --build ${build} --o "results_no_filter"