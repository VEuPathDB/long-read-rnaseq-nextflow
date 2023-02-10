#!/usr/bin/env bash

set -euo pipefail

talon_abundance --db ${database}  -a ${annot_name} --build ${build} --o "results_no_filter"