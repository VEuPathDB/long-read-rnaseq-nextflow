#!/usr/bin/env bash

set -euo pipefail


sort -k1,1 -k4,4n ${gff} > ${params.build}_sorted.gff 
bgzip ${params.build}_sorted.gff 
tabix -p gff ${params.build}_sorted.gff.gz