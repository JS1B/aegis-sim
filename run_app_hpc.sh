#!/bin/bash
# HPC submission script for Aegis-Murmuration simulation
# Submit with: bsub < run_app_hpc.sh

#BSUB -J AegisMurmuration
#BSUB -q hpcintro
#BSUB -W 00:30
#BSUB -n 16
#BSUB -R "rusage[mem=2GB]"
#BSUB -R "span[hosts=1]"
#BSUB -o logs/%J.out
#BSUB -e logs/%J.err

# Load required modules
module load gcc/11.5.0

# Create output directories
mkdir -p output_data
mkdir -p logs

# Set OpenMP threads for future Stage 2 parallelization
export OMP_NUM_THREADS=$LSB_DJOB_NUMPROC

echo "=== HPC Job Info ==="
echo "Job ID: $LSB_JOBID"
echo "Host: $(hostname)"
echo "Cores: $LSB_DJOB_NUMPROC"
echo "OMP_NUM_THREADS: $OMP_NUM_THREADS"
echo ""

# Parameters: N_particles, timesteps, output_interval
# Increase particle count for HPC runs
N_PARTICLES=50000
TIMESTEPS=1000
OUTPUT_INTERVAL=10

# Execute the simulation
/bin/bash ./run_app.sh $N_PARTICLES $TIMESTEPS $OUTPUT_INTERVAL

echo ""
echo "=== Job Complete ==="

