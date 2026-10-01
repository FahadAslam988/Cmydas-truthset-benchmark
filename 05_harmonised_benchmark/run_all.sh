#!/usr/bin/env bash
source "${CMYDAS_CONFIG:?export CMYDAS_CONFIG=/path/to/config/paths.sh}"
# Runs steps 1-4 in order. Safe to re-run: finished hap.py outputs are skipped.
set -uo pipefail; cd "$(dirname "$0")"; source common.sh
log "===== Paper 1 harmonised benchmark START ====="
./s01_TS1_build_consensus_truth.sh && ./s02_TS1_happy_consensus.sh && ./s03_TS4_happy_simulation_allele.sh && ./s04_collect_results.sh
log "===== Paper 1 harmonised benchmark END ====="
