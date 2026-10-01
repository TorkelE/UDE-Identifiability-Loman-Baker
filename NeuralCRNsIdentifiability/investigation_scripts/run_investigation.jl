### Preparations ###

# Reads inputs

# Sets the correct folder and activates its environment.
begin
    cd("$(@__DIR__)/../..")
    using Pkg; Pkg.activate(".")
end

# Fetches packages.
begin
    using Distributed
    using NeuralCRNsIdentifiability
    using Optimization
    import Plots: savefig
    using Serialization
    using StableRNGs
end

# Reads the inputs. For now, we always want to print all information.
begin
    print = true
    run_tag, is_scan, scan_idxs, parallel_scan = interpret_input(ARGS)
end

# Print run information.
if print
    println("Initiates investigation run: $run_tag")
    isnothing(scan_idxs) || println("Scan index: $scan_idxs")
end

# Sets the correct folder and fetches run settings.
begin
    main_folder = "$(is_scan ? "Scans" : "SingleRuns")/$run_tag"
    run_folder = is_scan ? "$main_folder/scan_$scan_idxs" : main_folder
    settings = deserialize("$main_folder/settings_run.jls")
    is_scan && (settings = combined_settings(settings, deserialize("$run_folder/settings_override.jls")))
    data_file = "$(is_scan && isfile("$run_folder/measured_data.jls") ? run_folder : main_folder)/measured_data.jls"
    md = MeasuredData(data_file)
end

# Creates the random number generator.
rng = StableRNG(settings[:seed])

### Runs Investigation ###

# Fetches the CRN models.
ps_ranges, models = retrieve_models(settings; rng)


# Performs analysis of each of the model versions.
for rs in models
    print && println("\nStarting model analysis: $(getname(rs))")

    # Fits the model to the data.
    begin
        success_required = is_valid(rs) && !(haskey(settings, :chain))
        opt_prob = make_eoptprob(rs, md, settings[:n]; rng, final_prints = false, rerun_print = false, success_required)
        @time opt_sol = fit_model(opt_prob, settings[:n];  maxiters = settings[:maxiters])
        print && summarise_optimisation(opt_sol, md)
        opt_sol_serializable = convert(Vector{NeuralCRNsIdentifiability.ExpandedOptimizationSolution}, opt_sol.u)
        serialize("$run_folder/opt_sol_$(getname(rs)).jls", opt_sol_serializable)
        savefig(plot_all_fitted_simulation(opt_sol, md), "$run_folder/base_figures/fitted_sols_$(getname(rs)).png")
        for (i, f_deps) in enumerate(get_dependent_sps(rs))
            (length(f_deps) != 1) && continue
            file_ending = (length(get_dependent_sps(rs)) == 1) ? "" : "_$i"
            filename = "$run_folder/base_figures/fitted_funcs_$(getname(rs))$(file_ending).png"
            savefig(plot_fitted_funcs(rs, opt_sol, md; f_idx = i), filename)
        end
    end

    # Performs likelihood profile analysis.
    begin
        all_pl_sols = []
        for method in settings[:profiling_methods]
            @time pl_sols = compute_lprofiles(method, rs, opt_prob, opt_sol, md)
            push!(all_pl_sols, pl_sols)
        end
        serialize("$run_folder/pl_sols_$(getname(rs)).jls", all_pl_sols)
        if length(all_pl_sols) == 1
            savefig(plot_lprofiles(all_pl_sols[1], md; info_mode = true), "$run_folder/base_figures/likelihood_profs_$(getname(rs)).png")
        elseif length(all_pl_sols) > 1
            for (idx, pl_sol) in enumerate(all_pl_sols)
                savefig(plot_lprofiles(pl_sol, md; info_mode = true), "$run_folder/base_figures/likelihood_profs_$(getname(rs))_$idx.png")
            end
        end
    end
end
