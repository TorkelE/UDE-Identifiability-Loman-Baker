#!/bin/bash
#SBATCH --job-name=CRNUDE_scan_extended_sal_basic_scan
#SBATCH --output=Scans/extended_sal_basic_scan/scan_%a/output_%A.log
#SBATCH --error=Scans/extended_sal_basic_scan/scan_%a/error_%A.log
#SBATCH --array=1-9
#SBATCH --nodes=1
#SBATCH --partition=medium
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --ntasks-per-node=1

# Print this sub-job's task ID
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

# Run the scan.
module purge
module load Julia/1.11.3-linux-x86_64
julia NeuralCRNsIdentifiability/investigation_scripts/run_investigation.jl extended_sal_basic_scan $SLURM_ARRAY_TASK_ID parallel_scan
