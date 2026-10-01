# Read inputs.
begin
    run_tag = (length(ARGS) == 0) ? "extended_sal_tmp_tmp_scan" : only(ARGS)
    run_folder = "Scans/$run_tag"
    isdir(run_folder) || error("There seems to be no scan associated with the run tag $run_tag")
end
# Fetches required packages. Sets the correct environment.
begin
    cd("$(@__DIR__)/../..")
    using Pkg; Pkg.activate(".")
    using NeuralCRNsIdentifiability
end

# Finds the missing scans.
missing_scans = check_scan(run_tag; print_info = false)

# Clears the folders of the scans which needs to be re-run.
for scan_idx in missing_scans
    for file in readdir("$run_folder/scan_$scan_idx")
        if (file == "measured_data.jls") || (file == "settings_override.jls")
            continue
        elseif file == "base_figures"
            for fig in readdir("$run_folder/scan_$scan_idx/base_figures")
                (fig == "measured_data.png") && continue
                rm("$run_folder/scan_$scan_idx/base_figures/$fig")
            end
        else
            rm("$run_folder/scan_$scan_idx/$file")
        end
    end
end

# Updates the compute_run.sh file.
begin
    idx_str = ""
    for idx in missing_scans
        global idx_str
        idx_str = "$idx_str $idx"
    end
    rm("NeuralCRNsIdentifiability/investigation_scripts/compute_run.sh", force = true)
    open("NeuralCRNsIdentifiability/investigation_scripts/compute_run.sh", "a") do io
        file_content = """
        #!/usr/bin/env bash

        set -euo pipefail

        M=25                 # max concurrent jobs

        # Enforce "1 process, 2 threads" per job
        export JULIA_NUM_THREADS=2
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        export MKL_NUM_THREADS=1
        export BLIS_NUM_THREADS=1
        export NUMEXPR_NUM_THREADS=1

        echo "Starting parallel scan for: $run_tag"
        parallel -j "\$M" --bar --line-buffer '
            idx="{}"
            dir="Scans/$run_tag/scan_\${idx}"
            echo "Starting scan run: \${idx}" | tee -a "\$dir/output.log"
            julia "NeuralCRNsIdentifiability/investigation_scripts/run_investigation.jl" $run_tag "\${idx}" >> "\$dir/output.log" 2>> "\$dir/error.log" parallel_scan
            echo "Finished scan run: \${idx}" | tee -a "\$dir/output.log"
        ' :::$idx_str
        echo "Finished parallel scan for: $run_tag"
        """
        write(io, file_content)
        chmod("NeuralCRNsIdentifiability/investigation_scripts/compute_run.sh", 0o755)
    end
end

# Updates the hpc_run.sh file.
begin
    idx_str_hpc = ""
    for idx in missing_scans
        global idx_str_hpc
        idx_str_hpc = "$idx_str_hpc,$idx"
    end
    rm("NeuralCRNsIdentifiability/investigation_scripts/hpc_run.sh", force = true)
    open("NeuralCRNsIdentifiability/investigation_scripts/hpc_run.sh", "a") do io
        file_content = """
    #!/bin/bash
    #SBATCH --job-name=CRNUDE_scan_$run_tag
    #SBATCH --output=Scans/$run_tag/scan_%a/output_%A.log
    #SBATCH --error=Scans/$run_tag/scan_%a/error_%A.log
    #SBATCH --array=$idx_str_hpc
    #SBATCH --nodes=1
    #SBATCH --partition=medium
    #SBATCH --time=48:00:00
    #SBATCH --mem-per-cpu=8G
    #SBATCH --ntasks-per-node=1

    # Print this sub-job's task ID
    echo "My SLURM_ARRAY_TASK_ID: " \$SLURM_ARRAY_TASK_ID

    # Run the scan.
    module purge
    module load Julia/1.11.3-linux-x86_64
    julia NeuralCRNsIdentifiability/investigation_scripts/run_investigation.jl $run_tag \$SLURM_ARRAY_TASK_ID parallel_scan
    """
        write(io, file_content)
    end
end
