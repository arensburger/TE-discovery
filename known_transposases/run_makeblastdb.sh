#!/bin/bash
# bash ./run_makeblastdb.sh
set -euo pipefail

makeblastdb -in tnps_sequences.fa -dbtype prot
