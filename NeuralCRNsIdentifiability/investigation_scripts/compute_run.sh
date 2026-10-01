#!/usr/bin/env bash
#
# Runs all points of a scan on the local machine.
# Run from the repository root, after setup_investigation.jl has prepared the scan.

set -euo pipefail

# ---- settings ---------------------------------------------------------------
TAG="extended_sal_basic_scan"     # the scan to run
M="${M:-}"                        # max concurrent jobs (default: half the cores)
# -----------------------------------------------------------------------------

# Number of scan points, taken from the folders setup_investigation.jl created.
N=$(find "Scans/${TAG}" -maxdepth 1 -type d -name 'scan_*' 2>/dev/null | wc -l | tr -d ' ')
if [ "$N" -eq 0 ]; then
    echo "No scan points found in Scans/${TAG}." >&2
    echo "Run setup_investigation.jl for this tag first." >&2
    exit 1
fi

# Each job uses 2 Julia threads, so default to half the available cores.
if [ -z "$M" ]; then
    CORES=$( (command -v nproc >/dev/null && nproc) || sysctl -n hw.ncpu 2>/dev/null || echo 2 )
    M=$(( CORES / 2 )); [ "$M" -lt 1 ] && M=1
fi

# Enforce "1 process, 2 threads" per job.
export JULIA_NUM_THREADS=2
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "Starting scan: ${TAG} (${N} points, ${M} at a time)"
seq 1 "$N" | xargs -P "$M" -I{} bash -c '
    idx="$1"; tag="$2"
    dir="Scans/${tag}/scan_${idx}"
    echo "Starting scan run: ${idx}" | tee -a "$dir/output.log"
    julia "NeuralCRNsIdentifiability/investigation_scripts/run_investigation.jl" \
        "$tag" "$idx" parallel_scan >> "$dir/output.log" 2>> "$dir/error.log"
    echo "Finished scan run: ${idx}" | tee -a "$dir/output.log"
' _ {} "$TAG"
echo "Finished scan: ${TAG}"
