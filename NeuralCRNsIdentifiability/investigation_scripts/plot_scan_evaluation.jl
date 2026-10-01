# Reads the inputs.
begin
    run_tag = ARGS[1]
    scan_type = if length(filter(s -> s ∉ ["skip_plot", "skip_vals"], ARGS)) == 2
        filter(s -> s ∉ ["skip_plot", "skip_vals"], ARGS)[2]
    elseif contains(run_tag, "nn_arch_scan")
        scan_type = "nn_arch_scan"
    elseif contains(run_tag, "data_variation_scan")
        scan_type = "data_variation_scan"
    else
        error("Could not determine scan type.")
    end
    make_plot = !in("skip_plot", ARGS)
    extract_vals = !in("skip_vals", ARGS)
end

# Print starting message.
println("Plotting the $scan_type-type evaluation for scan: $run_tag.")

# Sets the correct directory and import packages.
begin
    cd("$(@__DIR__)/../..")
    using Pkg; Pkg.activate(".")
    using NeuralCRNsIdentifiability
    using Serialization
    import Plots: savefig
end

# Loads the scan, checks it, and plots the result.
begin
    #Loading and checking.
    if !isempty(check_scan(run_tag; print_info = false))
        check_scan(run_tag; print_info = true)
    end
    @time scan = load_investigation(run_tag)

    # Plots the scan evaluation and saves it.
    plt = if scan_type == "nn_arch_scan"
        make_plot && @time plot_combined_nn_arch_eval(scan)
        make_plot && @time plot_full_data_quality_eval(scan)
        if extract_vals
            @time vals = extract_nn_architecture_grid_data(scan)
            serialize("Scans/$(run_tag)/extracted_vals.jls", vals)
        end
    else scan_type == "data_variation_scan"
        make_plot && @time plot_full_data_quality_eval(scan)
        if extract_vals
            @time vals = extract_datavar_scan_evaluation(scan)
            serialize("Scans/$(run_tag)/extracted_vals.jls", vals)
        end
    end
    make_plot && check_mkdir("Scans/$(run_tag)/base_figures")
    make_plot && savefig(plt, "Scans/$(run_tag)/base_figures/scan_evaluation.png")
end

# Print ending message.
println("Finished generating plot.")
