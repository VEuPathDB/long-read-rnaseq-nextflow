#!/usr/bin/env bash

set -euo pipefail

talon_initialize_database --f ${annotation} --a ${annot_name} --g ${build}  --o ${build}_"talon"