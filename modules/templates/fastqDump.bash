#!/usr/bin/env bash

set -euo pipefail

if [ -s $HOME/.ncbi/user-settings.mkfg ]
then
    fastq-dump ${sra}
else
    mkdir -p $HOME/.ncbi
    cp /usr/bin/user-settings.mkfg $HOME/.ncbi/
    fastq-dump ${sra}
fi