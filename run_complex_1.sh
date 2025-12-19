#!/bin/bash

#SBATCH --job-name=abfe_complex
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --array=0-2
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
# If you want GPUs for this job, add something like:
#SBATCH --gres=gpu:v100x:1

mkdir -p logs

# Activate conda (You may change method to activate conda here)
source ~/bin/myconda
conda activate rbfe

# Seed list (indexed by SLURM_ARRAY_TASK_ID)
SEEDS=(2025 2026 2027)
SEED=${SEEDS[$SLURM_ARRAY_TASK_ID]}

PDB="data/protein/4w53_fixed_keepMBN.pdb"
LIG_SDF="data/ligands/MBN.sdf"
RESNAME="MBN"
OUTDIR="results/complex_abs/MBN"
platform="auto" # or explicit "CPU", "CUDA", "OpenCL"

# Decide platform based on GPU availability

echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
echo "Using SEED=${SEED}"
echo "Detected CUDA_VISIBLE_DEVICES='${CUDA_VISIBLE_DEVICES:-none}'"
echo "Running on platform=${platform}"

python rbfe/run_complex_absolute.py \
  --pdb "$PDB" \
  --ligand_sdf "$LIG_SDF" \
  --ligand_resname "$RESNAME" \
  --seed "$SEED" \
  --platform "$platform" \
  --outdir "$OUTDIR"

