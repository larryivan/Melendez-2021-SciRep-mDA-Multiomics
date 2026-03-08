#!/usr/bin/env bash
set -euo pipefail

PIPELINE_PATH="${RNASEQ_PIPELINE_PATH:-/workspace/nextflow/rnaseq}"

show_help() {
    cat <<'EOF'
RNA-seq pipeline container

Usage:
  docker run ... melendez-rnaseq:latest help
  docker run ... melendez-rnaseq:latest run [pipeline args]
  docker run ... melendez-rnaseq:latest config [extra nextflow config args]
  docker run ... melendez-rnaseq:latest nextflow [nextflow args]
  docker run ... melendez-rnaseq:latest bash

Behavior:
  - `run` executes: nextflow run /workspace/nextflow/rnaseq -profile standard ...
  - unknown leading arguments are treated as pipeline arguments for `run`
  - mount the repository at /workspace for normal use

Examples:
  docker run --rm --platform linux/amd64 -v "$PWD":/workspace -w /workspace \
    melendez-rnaseq:latest run --samplesheet nextflow/rnaseq/assets/test_one_sample.tsv \
    --outdir results/rnaseq_nf_test -resume

  docker run --rm --platform linux/amd64 -v "$PWD":/workspace -w /workspace \
    melendez-rnaseq:latest config
EOF
}

cmd="${1:-help}"

case "${cmd}" in
    help|-h|--help)
        show_help
        ;;
    run)
        shift || true
        exec nextflow run "${PIPELINE_PATH}" -profile standard "$@"
        ;;
    config)
        shift || true
        exec nextflow config "${PIPELINE_PATH}" -profile standard "$@"
        ;;
    nextflow)
        shift || true
        exec nextflow "$@"
        ;;
    bash|sh)
        exec "${cmd}"
        ;;
    *)
        exec nextflow run "${PIPELINE_PATH}" -profile standard "$@"
        ;;
esac
