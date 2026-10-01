### Investigation Loading ###

# A structure storing the outputs from a single run.
struct SavedRun
    run_tag::String
    settings::Dict{Symbol,Any}
    ps_ranges::Vector
    models::Tuple
    md::MeasuredData
    opt_sols::Vector
    pl_sols::Vector

    # Base constructor.
    function SavedRun(run_tag::String, settings::Dict{Symbol,Any}, ps_ranges::Vector,
            models::Tuple, md::MeasuredData, opt_sols::Vector, pl_sols::Vector)
        new(run_tag, settings, ps_ranges, models, md, opt_sols, pl_sols)
    end

    # Specialised constructor used for when loading a single run.
    function SavedRun(run_tag::String, run_folder::String)
        # Finds the settings, parameter ranges, models, and measured data.
        settings = deserialize(run_folder *"/settings_run.jls")
        ps_ranges, models = retrieve_models(settings; rng = StableRNG(settings[:seed]))
        md = MeasuredData("SingleRuns/$run_tag/measured_data.jls")

        # Returns the output structure.
        SavedRun(run_tag, run_folder, settings, ps_ranges, models, md)
    end

    # Specialised constructor used for when loading a run which is part of a scan.
    function SavedRun(run_tag::String, run_folder::String, settings::Dict{Symbol,Any},
            ps_ranges::Vector, models::Tuple, md::MeasuredData)
        opt_sols = [deserialize_nothing("$run_folder/opt_sol_$name.jls") for name in getname.(models)]
        pl_sols =[deserialize_nothing("$run_folder/pl_sols_$name.jls") for name in getname.(models)]

        # Creates the correct data.
        new(run_tag, settings, ps_ranges, models, md, opt_sols, pl_sols)
    end
end

# Tries to load a file, but if it doesn't exist, returns `nothing` instead of an error.
function deserialize_nothing(filename::String)
    return isfile(filename) ? deserialize(filename) : nothing
end

# A structure storing the outputs from a scan.
struct SavedScan
    run_tag::String
    main_folder::String
    scan_settings::Vector
    runs::Vector{SavedRun}

    function SavedScan(run_tag::String, main_folder::String)
        # Finds all the scans.
        n_scans = length(filter(startswith("scan_"), readdir(main_folder)))

        # Loads base settings, and the parameter ranges and models.
        # Finds the settings, parameter ranges, models, and measured data.
        settings_general = deserialize(main_folder *"/settings_run.jls")
        ps_ranges, models = retrieve_models(settings_general; rng = StableRNG(settings_general[:seed]))

        # Load all the measured data (can be a single one or multiple).
        md = if isfile("$main_folder/measured_data.jls")
            fill(MeasuredData("$main_folder/measured_data.jls"), n_scans)
        else
            [MeasuredData("$main_folder/scan_$scan_idx/measured_data.jls") for scan_idx in 1:n_scans]
        end

        # Creates the unique settings for each scan.
        settings = [
            combined_settings(settings_general, deserialize("$main_folder/scan_$scan_idx/settings_override.jls"))
            for scan_idx in 1:n_scans
        ]

        # Finds all the runs.
        runs = [SavedRun(run_tag, main_folder * "/scan_$scan_idx", settings[scan_idx], ps_ranges, models, md[scan_idx])
            for scan_idx in 1:n_scans
        ]

        # Creates the vector of all scan settings.
        scan_settings = deserialize("$main_folder/settings_scan.jls")

        # Returns the full structure.
        new(run_tag, main_folder, scan_settings, runs)
    end
end

# For a run tag, loads the corresponding saved run or scan.
function load_investigation(run_tag::String)
    # Find the correct folder.
    (sum([run_tag in readdir("Scans"), run_tag in readdir("SingleRuns")]) != 1) && error("Run tag '$run_tag' must exist exactly once in either the SingRuns or Scans folder.")
    is_scan = (run_tag in readdir("Scans"))
    main_folder = "$(is_scan ? "Scans" : "SingleRuns")/$run_tag"

    # Returns the correct output depending on whether we have a scan or not.
    is_scan ? SavedScan(run_tag, main_folder) : SavedRun(run_tag, main_folder)
end

### Normal Run Plotting ###

# Creates a summaries plots of a single run.
function make_summarisation_plot(sr::SavedRun; x_pixs = 450, y_pixs = 350)
    # Creates all plots.
    fitting_plots = [plot_all_fitted_simulation(sr, idx)  for idx in eachindex(sr.models)]
    fitted_funcs_plots = [plot_fitted_funcs(sr, idx) for idx in eachindex(sr.models)]
    lp_prof_plots = [plot_lprofiles(sr, idx) for idx in eachindex(sr.models)]

    # Finds the grid size.
    n = sum(length.([fitting_plots, fitted_funcs_plots, lp_prof_plots]))
    n1 = argmin(2*n1 + 2*ceil(Int64, n/n1) for n1 in 1:ceil(Int64, sqrt(n)))
    n2 = ceil(Int64, n / n1)

    # Makes the final plot.
    plot(fitting_plots..., fitted_funcs_plots..., lp_prof_plots...; layout = (n1,n2), size = (n1*x_pixs, n2*y_pixs))
end

# Plots all the likelihood profiles in a single plot
function plot_combined_lprofiles(sr::SavedRun, p; model_idxs = eachindex(sr.models),
        colors = model_type_color.(sr.models), lws = fill(6, length(sr.models[model_idxs])),
        las = fill(0.7, length(sr.models[model_idxs])), kwargs...)
    plot()
    for (idx, model, color) in zip(model_idxs, sr.models[model_idxs], colors)
        plot_lprofile!(sr, idx, p; color, label = String(nameof(model)), lw = lws[idx], la = las[idx],
            plot_param = idx == model_idxs[end], plot_conf = idx == model_idxs[end], kwargs...)
    end
    plot!()
end

# Plots the combined likelihood profiles for all parameters.
function plot_all_combined_lprofiles(sr::SavedRun;  model_idxs = eachindex(sr.models), x_pixs = 450,
        y_pixs = 350, kwargs...)
    plots = [plot_combined_lprofiles(sr, p; model_idxs, kwargs...) for p in  MTK.getname.(parameters(sr.md.rs))]

    # Finds the grid size.
    n = length(plots)
    n1 = argmin(2*n1 + 2*ceil(Int64, n/n1) for n1 in 1:ceil(Int64, sqrt(n)))
    n2 = ceil(Int64, n / n1)

    # Makes the final plot.
    plot(plots...; layout = (n1,n2), size = (n1*x_pixs, n2*y_pixs))
end

# Plots the functional rates for a parameter and a nn model.
function plot_combined_fitted_funcs(sr::SavedRun; lw_params = 0.7, lw_nn = 0.2, kwargs...)
    param_idx = findfirst(contains(String(MTK.getname(model)), "params") for model in sr.models)
    nn_idx = findfirst(contains(String(MTK.getname(model)), "nn") for model in sr.models)
    plot_fitted_funcs(sr, param_idx; plot_true_rate = false, func_label = "Hill rate", la_bounds = lw_params, kwargs...)
    plot_fitted_funcs!(sr, nn_idx; func_label = "NN rate", la_bounds = lw_nn, kwargs...)
end

# Plots the fitted solution for each of the model types.
function plot_fitted_sims_combined(sr::SavedRun; cond = 1, label_data = true,
        show_nonmeasured = true, label_nonmeasured = true,
        lw_true = 9, la_true = 0.7, alpha_true = 0.6, ms_true = 4,
        colors = plotting_colors_alt,
        plot_models = eachindex(sr.models), linestyles = [:dash, :dot, :dashdot], lws = fill(5, length(sr.models)),
        las = fill(0.6, length(sr.models)), kwargs...)
    md = sr.md
    if num_measured(md) > length(colors)
        error("The number of measured species ($(num_measured(md))) exceeds the number of available plotting colors ($(length(plotting_colors))).")
    end
    plot()

    # If we wish to plot the non-measured species, plot their trajectories (faintly).
    if show_nonmeasured
        for (i, m) in enumerate(get_non_measured(md))
            color = colors[findfirst(isequal(m), MTK.getname.(species(md.rs)))][1]
            label = label_nonmeasured ? "$m (true)" : ""
            plot!(md.sols_true[cond]; idxs = m, color, label, lw = lw_true, la = la_true)
        end
    end

    # Plots the measured trajectories and the corresponding measurements.
    for (i, m) in enumerate(md.measured_sps)
        color1 = colors[findfirst(isequal(m), MTK.getname.(species(md.rs)))][1]
        color2 = colors[findfirst(isequal(m), MTK.getname.(species(md.rs)))][2]
        plot!(md.sols_true[cond]; idxs = m, color = color1, label = "$m (true)", lw = lw_true, la = la_true)
        label = label_data ? "$m (measured)" : ""
        plot!(md.measurement_ts, md.measurement_us[cond][m]; seriestype = :scatter,
            color = color2, label, alpha = alpha_true, ms = ms_true)
    end

    for (idx, model_idx) in enumerate(plot_models)
        model = sr.models[model_idx]
        lw = lws[idx]; la = las[idx]
        color_idx = idx + 2
        linestyle = linestyles[idx]
        # begin
        #     modelname = string(MTK.nameof(model))
        #     color_idx = if contains(modelname, "known")
        #         3
        #     elseif contains(modelname, "params")
        #         4
        #     elseif contains(modelname, "nn")
        #         5
        #     else
        #         Inf
        #     end
        # end
        opt_sol = sr.opt_sols[idx]

        # Simulates the fitted parameter set. If requested, plots the non-measured species.
        fitted_osol = get_fitted_osol(opt_sol; cond)
        if show_nonmeasured
            for (i, m) in enumerate(get_non_measured(md))
                color = colors[findfirst(isequal(m), MTK.getname.(species(md.rs)))][color_idx]
                label = label_nonmeasured ? "$m (fitted)" : ""
                plot!(fitted_osol; linestyle, idxs = m, color, label, lw, la)
            end
        end

        # Loops through each measured species and plots its fit. Returns the plot.
        for (i, m) in enumerate(md.measured_sps)
            color = colors[findfirst(isequal(m), MTK.getname.(species(md.rs)))][color_idx]
            plot!(fitted_osol; idxs = m, color, linestyle, label = "$m (fitted)", lw, la)
        end
    end

    # Performs final adjustments and returns the plot.
    plot!(; xlimit = (md.measurement_ts[1], md.measurement_ts[end]), ylimit = (0.0, Inf), kwargs...)
end

# Plot some condensed information in a single plot.
function plot_condensed_info(sr::SavedRun; x_pixs = 450, y_pixs = 350, max_num = 50)
    plots = [plot_fitted_sims_combined(sr), plot_combined_fitted_funcs(sr; max_num)]
    foreach(p -> push!(plots, plot_all_combined_lprofiles(sr)), MTK.getname.(parameters(sr.md.rs)))

    # Finds the grid size.
    n = length(plots)
    n1 = argmin(2*n1 + 2*ceil(Int64, n/n1) for n1 in 1:ceil(Int64, sqrt(n)))
    n2 = ceil(Int64, n / n1)

    # Makes the final plot.
    plot(plots...; layout = (n1,n2), size = (n2*x_pixs, n1*y_pixs),
        bottom_margin = 7Plots.mm, left_margin = 4Plots.mm)
end

# Helper function which, for a type of color, finds a corresponding color.
function model_type_color(model)
    modelname = string(MTK.nameof(model))
    if contains(modelname, "known")
        return :darkorange
    elseif contains(modelname, "nn_params") || contains(modelname, "nn_mondec_bounded") || contains(modelname, "nn_moninc_bounded")
        return :hotpink2
    elseif contains(modelname, "params_nn") || contains(modelname, "nn_mondec") || contains(modelname, "nn_moninc")
        return :turquoise
    elseif contains(modelname, "params")
        return :purple
    elseif contains(modelname, "nn")
        return :green
    elseif contains(modelname, "wrong")
        return :cyan
    else
        return :black
    end
end

### Specialised Scan Plots: NN Architecture Evaluation ###

# Plots the combined neural network architecture evaluation (effect of size on fit and overfitting).
function plot_combined_nn_arch_eval(ss::SavedScan)
    fit_grid = plot_nn_arch_fit_eval(ss; title = "Goodness-of-fit")
    overfit_grid = plot_nn_arch_overfit_eval(ss; title = "Degree of overfitting")
    return plot(fit_grid, overfit_grid; layout = (2,1), size = (600, 900), left_margin = 12Plots.mm, bottom_margin = 12Plots.mm)
end

# Plots the evaluation of how neural network size affects the fit quality.
function plot_nn_arch_fit_eval(ss::SavedScan; kwargs...)
    fitness_grid = generate_fitness_grid(ss)
    lengs, widths = get_nn_size_grid(ss)
    heatmap(lengs, widths, fitness_grid; xguide = "Num intermediary layers", yguide = "Layer width", kwargs...)
end

# Generates the fitness grid.
function generate_fitness_grid(ss::SavedScan)
    lengs, widths = get_nn_size_grid(ss)
    fitness_grid = [
        objective_nothing(ss.runs[nn_grid_idx_to_scan_idx(leng, width, ss)].opt_sols[end])
        for width in widths, leng in lengs
    ]

    # Normalises the fitness grid by the true fitness.
    fitness_grid ./= true_loss(ss.runs[1].md)
    return fitness_grid
end

# Special objective which returns NaN if the input is nothing.
objective_nothing(osol) = isnothing(osol) ? NaN : objective(osol)

# Plots the evaluation of how neural network size affects overfitting.
function plot_nn_arch_overfit_eval(ss::SavedScan; kwargs...)
    overfitting_grid = generate_overfitting_grid(ss)
    lengs, widths = get_nn_size_grid(ss)
    heatmap(lengs, widths, overfitting_grid; xguide = "Num intermediary layers", yguide = "Layer width", kwargs...)
end

# Generates the overfitting grid.
function generate_overfitting_grid(ss::SavedScan)
    # Computes the overfitting grid across all neural network architectures.
    lengs, widths = get_nn_size_grid(ss)
    return [
        compute_overfitting(ss.runs[nn_grid_idx_to_scan_idx(leng, width, ss)])
        for width in widths, leng in lengs
    ]
end

# Retrieves the grid of neural network sizes that a scan covers.
function get_nn_size_grid(ss::SavedScan)
    nn_sizes = [get_nn_size(setting[:chain]) for setting in ss.scan_settings]
    min_leng, max_leng = extrema(first.(nn_sizes))
    min_width, max_width = extrema(last.(nn_sizes))
    return (min_leng:max_leng, min_width:max_width)
end

# For a single neural network, finds its size.
function get_nn_size(nn)
    leng = length(nn.layers) - 2
    width = nn.layers[1].out_dims
    return (leng, width)
end

# From a grid index in the neural network architecture, get the flattened scan idx.
function nn_grid_idx_to_scan_idx(leng::Int64, width::Int64, ss::SavedScan)
    return findfirst(
        (length(set[:chain].layers) == leng + 2) && (set[:chain].layers[1].out_dims == width)
        for set in ss.scan_settings
    )
end

# Computes a measure of overfitting for a run (by comparing the found solution to the true solution).
# Can likely improve measure of sol difference (e.g. by using proper integral computation tool).
function compute_overfitting(sr::SavedRun, model_idx::Int64 = length(sr.models))
    # If the optimisation solution was not computed, returns a NaN.
    isnothing(sr.opt_sols[model_idx]) && return NaN

    # If the fitted run have worse objective value, there is no overfitting.
    (objective(sr.opt_sols[model_idx]) > true_loss(sr.md)) && return 0.0

    # Computes the total different between the two functions.
    overfitting = 0
    for cond in 1:numconditions(sr.md)
        fitted_osol = get_fitted_osol(sr.opt_sols[model_idx]; cond)
        true_osol = sr.md.sols_true[cond]
        t_vals = range(sr.md.measurement_ts[1], sr.md.measurement_ts[end], length = 100)
        for sp in sr.md.measured_sps
            overfitting += sum(abs, true_osol(t_vals; idxs = sp) .- fitted_osol(t_vals; idxs = sp))
        end
    end
    return overfitting
end

### Data Quality Scan Evaluation ###

# Plots the full evaluation of the effect of the data for all models and types of data.
function plot_full_data_quality_eval(ss::SavedScan; model_idxs = eachindex(ss.runs[1].models),
        x_pixs = 770, y_pixs = 360, left_margin = 12Plots.mm, bottom_margin = 12Plots.mm, kwargs...)
    plots = [plot_mean_func_diff(ss; model_idxs, kwargs...)]
    isempty(ss.runs[1].pl_sols) || push!(plots, plot_profile_likelihood_intervals(ss; model_idxs, kwargs...))
    if length(ss.runs[1].md.measured_sps) < length(species(ss.runs[1].md.rs))
        push!(plots, plot_nonmeasured_mean_err(ss; model_idxs, kwargs...))
    end
    return plot(plots...; layout = (length(plots),1), size = (x_pixs, length(plots)*y_pixs),
        left_margin, bottom_margin)
end

# Plots the profile likelihood confidence interval trend.
function plot_profile_likelihood_intervals(ss::SavedScan; model_idxs = eachindex(ss.runs[1].models), kwargs...)
    return plot_fit_property_of_data_quality(get_profile_likelihood_intervals, ss; model_idxs, yguide = "Mean conf, int.", kwargs...)
end

# Plots the non-measured mean species error trend.
function plot_nonmeasured_mean_err(ss::SavedScan; model_idxs = eachindex(ss.runs[1].models), kwargs...)
    return plot_fit_property_of_data_quality(get_nonmeasured_mean_err, ss; model_idxs, yguide = "Non-measured mean error", kwargs...)
end

# Plots the mean fitted function difference trend.
function plot_mean_func_diff(ss::SavedScan; model_idxs = eachindex(ss.runs[1].models), kwargs...)
    return plot_fit_property_of_data_quality(get_mean_func_diff, ss; model_idxs, yguide = "Mean function difference", kwargs...)
end

# For a selected property, performs its plot for all models (in a single plot).
function plot_fit_property_of_data_quality(eval_func, ss::SavedScan; model_idxs = eachindex(ss.runs[1].models), kwargs...)
    plot()
    foreach(model_idx -> plot_fit_property_of_data_quality!(eval_func, ss, model_idx; kwargs...), model_idxs)
    plot!()
end

# For a selected property, performs its plot for a single model.
plot_fit_property_of_data_quality(args...; kwargs...) = (plot(); plot_fit_property_of_data_quality!(args...; kwargs...);)
function plot_fit_property_of_data_quality!(eval_func, ss::SavedScan, model_idx::Int64;
        label = String(MTK.getname(ss.runs[1].models[model_idx])), color = model_type_color(ss.runs[1].models[model_idx]), kwargs...)
    scan_vals = get_scan_values(eval_func, ss::SavedScan, model_idx::Int64)
    plot!(eachindex(scan_vals), mean.(scan_vals); ribbon = std.(scan_vals), xguide = "Data quality", label, color, kwargs...)
end

# Gets the scan values for a specific scan and evaluation function.
function get_scan_values(eval_func, ss::SavedScan, model_idx::Int64)
    # Gets the scan values across the entire scan.
    leng, depth = get_scan_length_depth(ss)
    scan_vals = [
        [eval_func(ss.runs[m + (n-1)*depth], model_idx) for m in 1:depth]
        for n in 1:leng
    ]

    # Filters away any NaN values.
    scan_vals = [filter(!isnan, vals) for vals in scan_vals]
end

# Returns the data quality scan's length and depth.
function get_scan_length_depth(ss::SavedScan)
    leng = length(unique(delete!.(deepcopy(ss.scan_settings), :seed)))
    depth = Int64(length(ss.scan_settings)/leng)
    return (leng, depth)
end

# Retrieves the profile likelihood-based confidence interval for a parameter.
function get_profile_likelihood_intervals(sr::SavedRun, model_idx::Int64; p = fitted_parameters(sr.md.rs)[1])
    # If the profile was not computed, returns a NaN.
    isnothing(sr.pl_sols[model_idx]) && return NaN

    # Computes the profile likelihood derived confidence interval.
    pl_sol = only(sr.pl_sols[model_idx])[p]
    ovals_scaled = minimum(pl_sol.obj) .- pl_sol.obj
    idx_start = findfirst(ovals_scaled .> -1.92)
    idx_end = findlast(ovals_scaled .> -1.92)
    x_start = if idx_start == 1
        pl_sol.x[idx_start]
    else
        x1 = pl_sol.x[idx_start - 1]
        x3 = pl_sol.x[idx_start]
        y1 = ovals_scaled[idx_start - 1]
        y2 = -1.92
        y3 = ovals_scaled[idx_start]
        x1 + (x3-x1)*(y2-y1)/(y3-y1)
    end
    x_end = if idx_end == length(ovals_scaled)
        pl_sol.x[idx_end]
    else
        x1 = pl_sol.x[idx_end]
        x3 = pl_sol.x[idx_end + 1]
        y1 = ovals_scaled[idx_end]
        y2 = -1.92
        y3 = ovals_scaled[idx_end + 1]
        x1 + (x3-x1)*(y2-y1)/(y3-y1)
    end
    return x_end - x_start
end

# Retrieves the mean difference between the predicted and true values for the non-measured species.
function get_nonmeasured_mean_err(sr::SavedRun, model_idx::Int64)
    # If the solution was not computed, returns a NaN. Else, computes the mean error.
    isnothing(sr.opt_sols[model_idx]) && return NaN
    return mean(get_nonmeasured_err(opt_sol, sr.md) for opt_sol in sr.opt_sols[model_idx])
end

# Retrieves the difference between the predicted and true values for the non-measured species.
function get_nonmeasured_err(opt_sol::ExpandedOptimizationSolution, md::MeasuredData)
    tot_diff = 0
    for cond in 1:numconditions(md)
        fitted_osol = get_fitted_osol(opt_sol; cond)
        true_osol = md.sols_true[cond]
        t_vals = range(md.measurement_ts[1], md.measurement_ts[end], length = 100)
        for sp in get_non_measured(md)
            tot_diff += sum(abs, true_osol(t_vals; idxs = sp) .- fitted_osol(t_vals; idxs = sp))
        end
    end
    return tot_diff
end

# Retrieves the mean fitted function diff for a single model.
# Currently does L1 distance, also try L2 and Linf.
function get_mean_func_diff(sr::SavedRun, model_idx::Int64; f_idx::Int64 = 1, dens::Int64 = 200, compression_factor = 1.05)
    # If the solution was not computed, returns a NaN.
    isnothing(sr.opt_sols[model_idx]) && return NaN

    # Computes the mean difference between all fitted functions.
    fs_evaluated = evaluate_all_funcs(sr, model_idx; f_idx, dens, compression_factor)
    mean_diff = 0
    for i = 1:length(fs_evaluated), j = (i+1:length(fs_evaluated))
        mean_diff += sum(abs, fs_evaluated[i] .- fs_evaluated[j])
    end
    return mean_diff / (length(fs_evaluated) * (length(fs_evaluated) - 1) / 2)
end

# For a run, returns all the evaluated functions as a vector of vectors.
function evaluate_all_funcs(sr::SavedRun, model_idx::Int64; f_idx::Int64 = 1, dens::Int64 = 200, compression_factor = 1.05)
    # Extracts some stuff.
    md = sr.md
    rs = sr.models[model_idx]
    opt_sols = sr.opt_sols[model_idx]

    # Finds the function's input species. Ensures that only good solutions are considered.
    sp_in = only(get_dependent_sps(rs)[f_idx])
    opt_sols = filter(os -> os.objective < true_loss(md), opt_sols)

    # Error check.
    (length(get_dependent_sps(rs)[f_idx]) != 1)  && error("Currently, `evaluate_all_funcs` only supports functions of a single variable. In the future, this would be possible though.")

    # Find the evaluation range. Typically, we limit ourselves to values naturally occurring ground truth simulations.
    f_input_vals = vcat([max.(sol[sp_in], 0.0) for sol in md.sols_true]...)
    f_min, f_max = extrema(f_input_vals)
    f_min = compression_factor*f_min; f_max = f_max/compression_factor
    f_eval_vals = range(f_min; stop = f_max, length = dens)

    f_idx_is_nn(rs, f_idx)

    symbolic_rate = get_rates_of_interest(rs)[f_idx]
    [MTK.wrap(MTK.substitute(symbolic_rate, Dict(Pair.(fitted_parameters(rs), best_u(opt_sol)))))
        for opt_sol in opt_sols]

    # Compute the fitted function's value for all fits (different approach for nns and non-nns).
    fs_evaluated = if f_idx_is_nn(rs, f_idx)
        oprob_base, _, _, _, set_ps, _, _, _, _, _, _ = get_opt_prob_params(md, rs)
        [eval_nn(f_eval_vals, rs, f_idx_to_nn_idx(rs, f_idx), opt_sol.u, md, oprob_base, set_ps) for opt_sol in opt_sols]
    else
        symbolic_rate = get_rates_of_interest(rs)[f_idx]
        num_vals = [MTK.wrap(MTK.substitute(symbolic_rate, Dict(Pair.(fitted_parameters(rs), best_u(opt_sol)))))
            for opt_sol in opt_sols]
        [
            [
                MTK.value(MTK.substitute(num_val, Dict(getproperty(rs, sp_in) => sp_val)))
                for sp_val in f_eval_vals
            ]
            for num_val in num_vals
        ]
    end
end


### Data Quality Scan Helper Function ###

# Plots the progression of the cond information.
function plot_condensed_info_progression(ss::SavedScan, run_idxs; x_pixs = 400, y_pixs = 200, max_num = 100)
    plots = []
    foreach(i -> push!(plots, plot_condensed_info_reversed(ss.runs[i]; max_num)), run_idxs)
    plot(plots..., layout = (1, length(plots)), size = (length(plots)*x_pixs, 3*y_pixs))
end

# Plot some condensed information in a single plot. Plots it reversed
function plot_condensed_info_reversed(sr::SavedRun; x_pixs = 400, y_pixs = 200, max_num = 100)
    plots = [plot_fitted_sims_combined(sr), plot_combined_fitted_funcs(sr; max_num)]
    foreach(p -> push!(plots, plot_all_combined_lprofiles(sr)), MTK.getname.(parameters(sr.md.rs)))

    # Finds the grid size.
    n = length(plots)
    n1 = argmin(2*n1 + 2*ceil(Int64, n/n1) for n1 in 1:ceil(Int64, sqrt(n)))
    n2 = ceil(Int64, n / n1)

    # Makes the final plot.
    plot(plots...; layout = (n2, n1), label = "", size = (n1*x_pixs, n2*y_pixs), xguide = "", yguide = "", xticks = [], yticks = [])
end

### Value Extraction Functions ###

# For a data variation scan, loads in all the values (which can then be processed to create plots).
function extract_datavar_scan_evaluation(saved_scan)
    all_vals = [[
        get_scan_values(get_profile_likelihood_intervals, saved_scan, model_idx),
        get_scan_values(get_nonmeasured_mean_err, saved_scan, model_idx),
        get_scan_values(get_mean_func_diff, saved_scan, model_idx),
    ] for model_idx in 1:length(saved_scan.runs[1].models)]
    mean_vals = [[mean.(vals) for vals in model_vals] for model_vals in all_vals]
    std_vals = [[std.(vals) for vals in model_vals] for model_vals in all_vals]
    return mean_vals, std_vals
end

# Extracts NN architecture grid data.
function extract_nn_architecture_grid_data(saved_scan)
    fitness_grid = generate_fitness_grid(saved_scan)
    overfitting_grid = generate_overfitting_grid(saved_scan)
    return fitness_grid, overfitting_grid
end
