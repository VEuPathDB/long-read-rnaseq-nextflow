#!/usr/bin/env bash

set -euo pipefail
touch ${annotationName}
talon  --f ${config} \
       --db ${database} \
       --build ${build} \
       --o results