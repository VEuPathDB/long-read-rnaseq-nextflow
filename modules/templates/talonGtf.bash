#!/usr/bin/env bash

set -euo pipefail

talon_create_GTF --db ${database} -a ${annot_name} --build ${build} --o results