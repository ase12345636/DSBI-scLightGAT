#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cache_dir="$script_dir/.cache"
python_executable=/Group16T/common/ccuc/miniconda3_gb10/envs/scLightGAT/bin/python3

mkdir -p "$cache_dir/matplotlib" "$cache_dir/numba"
export PYTHONDONTWRITEBYTECODE=1
export MPLCONFIGDIR="$cache_dir/matplotlib"
export NUMBA_CACHE_DIR="$cache_dir/numba"
export SCIPY_ARRAY_API=1

# Read the single resource setting without importing the model or numeric
# libraries, then apply it before the production Python process starts.
cpu_workers="$(
    "$python_executable" -c \
        'import runpy, sys; print(runpy.run_path(sys.argv[1])["CPU_WORKERS"])' \
        "$script_dir/scLightGAT/scLightGAT/hparam.py"
)"
if ! [[ "$cpu_workers" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: CPU_WORKERS in hparam.py must be a positive integer." >&2
    exit 2
fi
export OMP_NUM_THREADS="$cpu_workers"
export OPENBLAS_NUM_THREADS="$cpu_workers"
export MKL_NUM_THREADS="$cpu_workers"
export BLIS_NUM_THREADS="$cpu_workers"
export VECLIB_MAXIMUM_THREADS="$cpu_workers"
export NUMEXPR_NUM_THREADS="$cpu_workers"
export NUMBA_NUM_THREADS="$cpu_workers"

cd "$script_dir"
exec "$python_executable" scLightGAT/main.py "$@"
