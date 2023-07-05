#!/usr/bin/env bash

set -euo pipefail

touch ${annotOut}

talon_create_GTF --db ${database} -a ${annotationName} --build ${build} --o ${build}_results