#!/usr/bin/env bash
set -euo pipefail
if [[ $# -ne 1 ]]; then
    echo "Usage: $0 /shared/path/flagger-v1.2.0.sif" >&2
    exit 2
fi
if [[ -e "$1" ]]; then
    echo "Output already exists: $1" >&2
    exit 1
fi
# Apptainer imports upstream's OCI image into SIF; no Docker daemon is used.
apptainer pull "$1" docker://mobinasri/flagger:v1.2.0
apptainer exec "$1" bash -c 'command -v bam2cov && command -v hmm_flagger'
