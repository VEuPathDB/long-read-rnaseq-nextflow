#!/usr/bin/env bash

set -euo pipefail

talon_abundance --db ${database} --whitelist ${wishList}  -a ${annotationName} --build ${build} --o "results"