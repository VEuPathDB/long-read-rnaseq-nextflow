#!/usr/bin/env bash

set -euo pipefail

talon  --f ${config} --db ${database} --build ${build} --o results