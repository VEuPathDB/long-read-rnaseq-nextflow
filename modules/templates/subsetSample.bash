#!/usr/bin/env bash

set -euo pipefail


subset_by_sample.py ${unFilteredCounts} ${filteredCounts}
    