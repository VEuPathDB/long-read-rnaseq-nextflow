#!/usr/bin/env bash

set -euo pipefail
touch ${annot_name}
talon  --f ${config} --db ${database} --build ${build} --o results