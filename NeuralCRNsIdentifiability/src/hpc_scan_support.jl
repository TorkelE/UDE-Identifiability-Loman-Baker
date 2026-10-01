### Input Reading ####

function interpret_input(args::Vector{String})
    # Checks if it is a hpc run.
    parallel_scan = ("parallel_scan" in args)
    if parallel_scan
        args = filter(s -> s != "parallel_scan", args)
    end

    # Error checks.
    length(args) > 2 && error("The was an unexpected number of input arguments ($length(args)).")
    (!isempty(args)) && (sum(check_int.(args)) > 1) && error("There was more than one input argument that was an integer.")
    (length(args) == 2) && (sum(check_int.(args)) == 0) && error("If there are two input arguments, one of them must be an integer.")

    # Attempts to read the active run file (for now, we require it to exist).
    active_run = deserialize("NeuralCRNsIdentifiability/investigation_scripts/active_run.jls")

    # Reads the run tag.
    run_tag_arg_idx = findfirst(!check_int, args)
    run_tag = (isnothing(run_tag_arg_idx) ? active_run[:run_tag] : args[run_tag_arg_idx])

    # Figures out if this is a scan or not.
    scan_idx_arg_idx = findfirst(check_int, args)
    isdir("Scans/" * run_tag) && isdir("SingleRuns/" * run_tag) && error("The run tag $run_tag exists both as a scan and as a single run.")
    !isdir("Scans/" * run_tag) && !isdir("SingleRuns/" * run_tag) && error("The run tag $run_tag exists neither as a scan and as a single run.")
    is_scan = isdir("Scans/" * run_tag)

    # If this is a scan, figures out the scan_idx.
    scan_idx = if is_scan
        isnothing(scan_idx_arg_idx) ? active_run[:scan_idxs] : parse(Int64, args[scan_idx_arg_idx])
    else
        nothing
    end

    # Under certain circumstances, updates the active run file.
    if isnothing(scan_idx_arg_idx) && is_scan
        tot_scans = active_run[:tot_scans]
        active_run[:scan_idxs] = (active_run[:scan_idxs] % tot_scans) + 1
        serialize("NeuralCRNsIdentifiability/investigation_scripts/active_run.jls", active_run)
    end

    # Returns the read input.
    return run_tag, is_scan, scan_idx, parallel_scan
end

# Checks if a string corresponds to an integer.
check_int(s::String) = !isnothing(tryparse(Int, s))
check_int(s::Number) = return isinteger(s)

### Scanning Functions ###

# Creates a dictionary containing the settings for a run.
function make_setting_dict(settings_file::String)
    # Runs the settings file. Creates the settings dictionary.
    include(settings_file)
    settings = Dict{Symbol, Any}()

    # Goes through all settings which *should* be declared in the file.
    settings[:num_procs] = num_procs
    settings[:ps_true] = ps_true
    settings[:retrieval_function] = retrieval_function
    settings[:u0s] = u0s
    settings[:tend] = tend
    settings[:sample_n] = sample_n
    settings[:measured_sps] = measured_sps
    settings[:n] = n
    settings[:σ_std] = σ_std
    settings[:maxiters] = maxiters
    settings[:profiling_methods] = profiling_methods

    # Goes through all settings which *might* be declared in the file (but might have defaults).
    settings[:seed] = (@isdefined seed) ? seed : 123
    settings[:noise_func] = (@isdefined noise_func) ? noise_func : scaled_normal_error

    # Returns the settings dictionary.
    return settings
end

# Creates a vector of dictionaries, each containing the override settings for a specific
# instance of the scan.
function make_scan_setting_dicts(settings_file::String)
    include(settings_file)
    return scan_iterator
end

# Performs all computations for a single run folder.
# Reads the investigation setting for that folder, and runs the investigation.
# An optional argument is used in case this is a sub-run of a large-scale scan.
function run_investigation(run_tag; scan_idx = nothing, re_generate_data = true, print = true,
        make_figures = false)

    # Depending on whether this is part of a scan or a single run, we use different settings and folders.
    if !isnothing(scan_idx)
        run_folder = "Scans/$run_tag/scan_$scan_idx"
        parent_folder = "Scans/$run_tag"
    else
        run_folder = "SingleRuns/$run_tag"
        parent_folder = run_folder
    end

    # Loads the run settings (and potential settings overrides).
    settings = deserialize("$parent_folder/investigation_settings.jls")::Dict{Symbol, Any}
    if isfile("$run_folder/investigation_settings_override.jls")
        override_settings = deserialize("$run_folder/investigation_settings_override.jls")
        foreach((k,v) -> settings[k] = v, override_settings)
    end

    # Sets a random number generator.
    rng = StableRNG(settings[:seed])

    # Retrieves the model.
    ps_ranges, rs_known, rs_params, rs_nn = settings[:retrieval_function](settings[:ps_true], rng)

    # Generates the data (or, potentially, loads it from a file).
    data_filepath = "$(scan_wide_data ? run_folder : parent_folder)/measured_data.jls"
    if re_generate_data
        md = MeasuredData(rs_known, ps_ranges, rng, settings)
        serialize(data_filepath, md)
    else
        md = MeasuredData(data_filepath)
    end

    # Loops through the different models, performing analysis for each.
    for (rs, rs_tag) in zip([rs_known, rs_params, rs_nn], ["known", "params", "nn"])
        print && println("\nStarting model analysis: $rs_tag")

        # Firs the model to the data.
        opt_prob = make_eoptprob(rs, md, settings[:n]; rng, final_prints = false)
        @time opt_sol = fit_model(opt_prob, settings[:n]; ealg = EnsembleThreads(), maxiters = settings[:maxiters])
        summarise_optimisation(opt_sol, md)
        opt_sol_serializable = convert(Vector{NeuralCRNsIdentifiability.ExpandedOptimizationSolution}, opt_sol.u)
        serialize("$run_folder/opt_sol_$rs_tag.jls", opt_sol_serializable)

        # Performs likelihood profile analysis.
        all_pl_sols = []
        for method in profiling_methods
            @time pl_sols = compute_lprofiles(method, rs, opt_prob, opt_sol, md)
            push!(all_pl_sols, pl_sols)
        end



        serialize("$run_folder/pl_sols_$rs_tag.jls", all_pl_sols)
    end

    # Creates the figures.
    if make_figures
        error("`run_investigation` does not yet support figure generation.")
    end

    # Prints end statement, analysis is now finished.
    print && println("Finished run: $run_tag")
end

# Generates a SLURM file for running a specific scan on the hpc.
# Also generate a .sh script for running things without SLURM on a compute machine.
function make_hpc_file(run_tag::String, n)
    # For slurm runs.
    rm("NeuralCRNsIdentifiability/investigation_scripts/hpc_run.sh", force = true)
    open("NeuralCRNsIdentifiability/investigation_scripts/hpc_run.sh", "a") do io
        file_content = """
    #!/bin/bash
    #SBATCH --job-name=CRNUDE_scan_$run_tag
    #SBATCH --output=Scans/$run_tag/scan_%a/output_%A.log
    #SBATCH --error=Scans/$run_tag/scan_%a/error_%A.log
    #SBATCH --array=1-$n
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

    # For compute runs.
    rm("NeuralCRNsIdentifiability/investigation_scripts/compute_run.sh", force = true)
    open("NeuralCRNsIdentifiability/investigation_scripts/compute_run.sh", "a") do io
        file_content = """
        #!/usr/bin/env bash

        set -euo pipefail

        N=$n                 # total jobs (1..N)
        M=25                 # max concurrent jobs

        # Enforce "1 process, 2 threads" per job
        export JULIA_NUM_THREADS=2
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        export MKL_NUM_THREADS=1
        export BLIS_NUM_THREADS=1
        export NUMEXPR_NUM_THREADS=1

        # Drive 1..N with at most M in parallel
        echo "Starting parallel scan for: $run_tag"
        seq 1 "\$N" | parallel -j "\$M" --bar --line-buffer '
            idx="{}"
            dir="Scans/$run_tag/scan_\${idx}"
            echo "Starting scan run: \${idx}" | tee -a "\$dir/output.log"
            julia "NeuralCRNsIdentifiability/investigation_scripts/run_investigation.jl" $run_tag "\${idx}" >> "\$dir/output.log" 2>> "\$dir/error.log" parallel_scan
            echo "Finished scan run: \${idx}" | tee -a "\$dir/output.log"
        '
        echo "Finished parallel scan for: $run_tag"
        """
        write(io, file_content)
        chmod("NeuralCRNsIdentifiability/investigation_scripts/compute_run.sh", 0o755)
    end
end


### Utility FUnctions ###

# Retrieves the CRN model set from a settings dictionary.
function retrieve_models(settings::Dict{Symbol, Any}; rng = Random.default_rng())
    if haskey(settings, :chain)
        return settings[:retrieval_function](settings[:ps_true], rng; chain = settings[:chain])
    else
        return settings[:retrieval_function](settings[:ps_true], rng)
    end
end

# If a directory doesn't exists, creates it.
check_mkdir(dir::String) = (isdir(dir) || mkdir(dir))

# Creates a merged dictionary (where the overwrite dictionary overwrites the base dictionary).
function combined_settings(base_settings, override_settings)
    settings_new = deepcopy(base_settings)
    foreach(entry -> settings_new[entry[1]] = entry[2], override_settings)
    return settings_new
end

### Scan Checking ###

# Evaluates to what extent a scan has been completed.
function check_scan(run_tag::String; print_info = true)
    # Finds the correct folder. Finds all of the relevant models.
    main_folder = "Scans/$run_tag"
    scan_length = length(filter(str -> startswith(str, "scan_"), readdir(main_folder)))
    settings = deserialize("$main_folder/settings_run.jls")
    _, models = retrieve_models(settings)
    model_names = [getname(rs) for rs in models]

    # Loops through all files that should have been created.
    incomplete_idxs = []
    for idx in 1:scan_length
        completed = true
        for model_name in model_names
            if !isfile("$main_folder/scan_$idx/opt_sol_$model_name.jls")
                print_info && println("Scan: $idx did not complete opt_sol_$model_name.jls")
                completed = false
            end
        end
        if !isempty(settings[:profiling_methods])
            for model_name in model_names
                if !isfile("$main_folder/scan_$idx/pl_sols_$model_name.jls")
                    print_info && println("Scan: $idx did not complete pl_sols_$model_name.jls")
                    completed = false
                end
            end
        end
        if !completed
            push!(incomplete_idxs, idx)
            print_info && println("\n")
        end
    end

    # Prints an output summary.
    print_info && print("Total incomplete scans: $(length(incomplete_idxs)) ($run_tag).")
    if !isempty(incomplete_idxs)
        print_info && println(" The following scans are incomplete:")
        print_info && print.(" ", incomplete_idxs)
        print_info && println()
    end
    return incomplete_idxs
end
