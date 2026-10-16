#!/bin/bash
#SBATCH --job-name=grand_1m
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=logs/sweep_%a.out
#SBATCH --error=logs/sweep_%a.err
#SBATCH --array=0-49%16

# ==============================================================================
# RhizAeon 1,000,000 Continuous Parameter Sweep Slurm Array Submission Script
# Total Trials: 50 tasks x 20,000 trials = 1,000,000 synthetic genomes
# Continuous distribution: K in {0, 1, 2, 3, 4} events across N in [4, 250],
# L in [300, 30,000] bp, divergence in [0.005, 0.35], HKY85 + discrete Gamma.
# ==============================================================================

TRIALS_PER_TASK=20000
START_TRIAL=$((SLURM_ARRAY_TASK_ID * TRIALS_PER_TASK))
OUT_DIR="results"
mkdir -p logs "${OUT_DIR}"
OUT_CSV="${OUT_DIR}/sweep_part_${SLURM_ARRAY_TASK_ID}.csv"

echo "================================================================================"
echo "Task ID:          ${SLURM_ARRAY_TASK_ID}"
echo "Node:             $(hostname)"
echo "CPUs Allocated:   ${SLURM_CPUS_PER_TASK}"
echo "Memory Allocated: 48G"
echo "Start Trial ID:   ${START_TRIAL}"
echo "Trials Count:     ${TRIALS_PER_TASK}"
echo "Output CSV:       ${OUT_CSV}"
echo "Start Time:       $(date)"
echo "================================================================================"

# Locate grand_sweep binary: PATH first, then a HyphAeon monorepo checkout.
# Build with: cargo run -p rhizaeon-tools --release --bin grand_sweep
SWEEP_BIN="$(which grand_sweep 2>/dev/null)"
if [ -z "${SWEEP_BIN}" ] && [ -n "${HYPHAEON_ROOT:-}" ]; then
    SWEEP_BIN="${HYPHAEON_ROOT}/rust/target/release/grand_sweep"
fi

if [ -z "${SWEEP_BIN}" ] || [ ! -x "${SWEEP_BIN}" ]; then
    echo "ERROR: grand_sweep binary not found."
    echo "  Build it in the HyphAeon monorepo: cargo build -p rhizaeon-tools --release"
    echo "  then put grand_sweep on PATH or set HYPHAEON_ROOT to the checkout."
    exit 1
fi

"${SWEEP_BIN}" \
    --start-trial ${START_TRIAL} \
    --trials ${TRIALS_PER_TASK} \
    --threads ${SLURM_CPUS_PER_TASK} \
    --output "${OUT_CSV}" \
    --seed $((42 + SLURM_ARRAY_TASK_ID * 101))

EXIT_CODE=$?
echo "================================================================================"
echo "End Time:         $(date)"
echo "Exit Code:        ${EXIT_CODE}"
echo "================================================================================"

exit ${EXIT_CODE}
