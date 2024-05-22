#!/usr/bin/env bash

set -euo pipefail

talon_initialize_database --f ${annotation} \
                          --a ${annotationName} \
                          --g ${build}  \
                          --o ${build}