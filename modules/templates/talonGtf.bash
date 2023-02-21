#!/usr/bin/env bash

set -euo pipefail

touch ${annotOut}

talon_create_GTF --db ${database} -a ${annot_name} --build ${build} --o ${build}_results