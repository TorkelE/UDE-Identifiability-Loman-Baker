# Functions for plotting the results of various types of analysis.

### Simulation Plotting ###

# Plots the data for all conditions.
function plot_all_data(md::MeasuredData)
    # Computes the best plot grid.
    n = numconditions(md)
    n1 = argmin(2*n1 + 2*ceil(Int64, n/n1) for n1 in 1:ceil(Int64, sqrt(n)))
    n2 = ceil(Int64, n / n1)

    # Plots the measured data (and the real simulations they are take from).
    plts = [plot_data(md; cond = i) for i in 1:numconditions(md)]
    plot(plts...; layout = (n1, n2), size = (n2*500, n1*350), bottom_margin = 10Plots.mm)
end

# Plots a simulated dataset.
plot_data(args...; kwargs...) = (plot(); plot_data!(args...; kwargs...);)
function plot_data!(md::MeasuredData, args...; cond = 1, show_nonmeasured = true,
        label_nonmeasured = show_nonmeasured, label_data = true, ms = 4, kwargs...)
    # Error checks.
    if num_measured(md) > length(plotting_colors)
        error("The number of measured species ($(num_measured(md))) exceeds the number of available plotting colors ($(length(plotting_colors))).")
    end

    # If we wish to plot the non-measured species, plot their trajectories (faintly).
    if show_nonmeasured
        for (i, m) in enumerate(get_non_measured(md))
            color = get_color(num_measured(md) + i)[1]
            label = label_nonmeasured ? "$m (true)" : ""
            plot!(md.sols_true[cond]; idxs = m, color, label, lw = 4, la = 0.5)
        end
    end

    # Plots the measured trajectories and the corresponding measurements.
    for (i, m) in enumerate(md.measured_sps)
        plot!(md.sols_true[cond]; idxs = m, color = plotting_colors[i][1], label = "$m (true)", lw = 7, la = 0.8)
        label = label_data ? "$m (measured)" : ""
        plot!(md.measurement_ts, md.measurement_us[cond][m]; seriestype = :scatter,
            color = plotting_colors[i][2], label = label, alpha = 0.6, ms)
    end

    # Performs final adjustments and returns the plot.
    plot!(; xlimit = (md.measurement_ts[1], md.measurement_ts[end]), ylimit = (0.0, Inf), kwargs...)
end

# Plots the fitted solution for all conditions.
function plot_all_fitted_simulation(opt_sol, md; label_nonmeasured = true, ms = 4, kwargs...)
    # Computes the best plot grid.
    n = numconditions(md)
    n1 = argmin(2*n1 + 2*ceil(Int64, n/n1) for n1 in 1:ceil(Int64, sqrt(n)))
    n2 = ceil(Int64, n / n1)

    # Plots the fits for all conditions.
    plts = [plot_fitted_simulation(opt_sol, md; cond = i, label_nonmeasured, ms, kwargs...) for i in 1:numconditions(md)]
    plot(plts...; layout = (n1, n2), size = (n2*500, n1*350), bottom_margin = 10Plots.mm)
end

# Dispatch for when a `SavedRun` object is provided.
function plot_all_fitted_simulation(sr::SavedRun, model_idx::Int64; kwargs...)
    return plot_all_fitted_simulation(sr.opt_sols[model_idx], sr.md; kwargs...)
end


# Plots a fitted simulation and compares it to the data.
plot_fitted_simulation(args...; kwargs...) = (plot(); plot_fitted_simulation!(args...; kwargs...);)
function plot_fitted_simulation!(opt_sol, md::MeasuredData, args...;
        cond = 1, show_nonmeasured = true, label_nonmeasured = true, ms = 4, kwargs...)
    # Plots the measured data.
    plot_data!(md; show_nonmeasured, label_nonmeasured = true, cond, ms)

    # Simulates the fitted parameter set. If requested, plots the non-measured species.
    fitted_osol = get_fitted_osol(opt_sol; cond)
    if show_nonmeasured
        for (i, m) in enumerate(get_non_measured(md))
            color = get_color(num_measured(md) + i)[3]
            label = label_nonmeasured ? "$m (fitted)" : ""
            plot!(fitted_osol; idxs = m, color, label, lw = 3, la = 1.0, linestyle = :dash)
        end
    end

    # Loops through each measured species and plots its fit. Returns the plot.
    for (i, m) in enumerate(md.measured_sps)
        plot!(fitted_osol; idxs = m, color = plotting_colors[i][3], linestyle = :dash, label = "$m (fitted)", lw = 3, la = 0.9)
    end
    plot!(; kwargs...)
end

### Identifiability Evaluation Plots ###

# Plots a single likelihood profile.
# Should come up with a nice way to automatically set the xguide to the parameter.
plot_lprofile(args...; kwargs...) = (plot(); plot_lprofile!(args...; kwargs...);)
function plot_lprofile!(pl_sol::LikelihoodProfile, md::MeasuredData; xguide = "", p = nothing,
        zoom = true, color = 1, la = 0.5, yguide = "loss value", label = "", yaxis = :identity,
        no_marker = false, info_mode = false, (yflip = yaxis == :log10), plot_conf = !yflip,
        plot_param = true, true_p_lw = 4, conf_lw = 3.0, ymin = -6.0, lw = 6, kwargs...)

    # Rescales the values. Determines finite ones and min and max values.
    ovals_scaled = minimum(pl_sol.obj) .- pl_sol.obj
    yflip && (ovals_scaled .= -ovals_scaled) + 1.0
    finites = isfinite.(ovals_scaled)
    (minval, maxval) = extrema(ovals_scaled[finites])

    # Handles potentially infinite points. Also handles markers for retcodes.
    if info_mode
        x = pl_sol.x
        y = [isfinite(o) ? o : (yflip ? maxval : minval) for o in ovals_scaled]
        markercolor = [isfin ? :orange2 : :red for isfin in finites]
        markershape = fill(:circle, length(x))
        ms = [7; [2.5 for x in pl_sol.x[2:end-1]]; 7]
        (pl_sol.retcodes[1] in [:Identifiable, :NonIdentifiable]) || (markershape[1] = :star4;)
        (pl_sol.retcodes[2] in [:Identifiable, :NonIdentifiable]) || (markershape[end] = :star4;)

        # Special marker for the MLE point.
        markershape[pl_sol.mle_idx] = :diamond
        ms[pl_sol.mle_idx] = 5
    else
        finites = isfinite.(pl_sol.obj)
        x = pl_sol.x[finites]
        y = ovals_scaled[finites]
        markercolor = fill(color, length(x))
        markershape = fill(:circle, length(x))
        ms = fill(2.5, length(x))
    end

    # Plots the profile.
    xlimit = (x[1] - abs(x[1])*0.02, x[end] + abs(x[end])*0.02)
    kwargs = merge_kwargs(kwargs; xlimit, yaxis, xguide, yguide, label, color, legend = :none, lw, la)
    no_marker || (kwargs = merge_kwargs(kwargs; ms))
    plot!(x, y; kwargs...)

    # If we know the designated parameter, we can do some additional stuff.
    if !isnothing(p) && plot_param
        xguide = MTK.getname(p)
        xlimit = md.ps_ranges[MTK.getname(p)]
        true_val = md.ps_true[MTK.getname(p)]
        plot!([true_val, true_val], collect(ylims(plot!())); xlimit, label = "$p (true)",
            color = :black, lw = true_p_lw, la = 0.5, linestyle = :dot, xguide = p)
    end

    # If requested, we mark the 95% confidence interval with a line.
    if plot_conf
        plot!(collect(xlimit), fill(-1.92, 2), color = :black, lw = conf_lw, la = 0.5,  linestyle = :dash, label = "95% CI")
    end

    # Likely we wish to zoom in on the relevant area. Returns plot.
    ylimit = (zoom ? (ymin, 0.2) : ylims())
    plot!(; ylimit, kwargs...)
end

# Dispatch for when a `SavedRun` object is provided.
function plot_lprofile!(sr::SavedRun, model_idx::Int64, p::Symbol; profile_method_idx = 1, kwargs...)
    pl_sol = sr.pl_sols[model_idx][profile_method_idx][sr.models[model_idx].var_to_name[p]]
    md = sr.md
    return plot_lprofile!(pl_sol, md; p, kwargs...)
end

# Plots all likelihood profiles. Makes sure to plot them in as nice a grid as possible.
 function plot_lprofiles(pl_sols::Dict{<:Any,LikelihoodProfile}, md::MeasuredData; core_ps_only = true, kwargs...)
    # Filters out parameters we do not want to plot (typically for the case of a parametrised function).
    plotted_ps = core_ps_only ? filter(p -> any(isequal(p), parameters(md.rs)), keys(pl_sols)) : keys(pl_sols)

    # Computes the best plot grid.
    n = length(plotted_ps)
    n1 = argmin(2*n1 + 2*ceil(Int64, n/n1) for n1 in 1:ceil(Int64, sqrt(n)))
    n2 = ceil(Int64, n / n1)

    # Plots all the profiles.
    plts = [plot_lprofile(pl_sols[p], md; yguide = "", p, kwargs...) for p in plotted_ps]
    plot(plts...; layout = (n1, n2), size = (n2*500, n1*350), bottom_margin = 10Plots.mm)
end

# Dispatch for when a `SavedRun` object is provided.
function plot_lprofiles(sr::SavedRun, model_idx::Int64; profile_method_idx = 1, kwargs...)
    return plot_lprofiles(sr.pl_sols[model_idx][profile_method_idx], sr.md; kwargs...)
end

### Other Plots ###

# Enables us to layer plots of fitted functions.
plot_fitted_func(args...; kwargs...) = (plot(); plot_fitted_func!(args...; kwargs...);)
plot_fitted_funcs(args...; kwargs...) = (plot(); plot_fitted_funcs!(args...; kwargs...);)

# For a non-ensemble optimization solution, plots the fitted functions (potential neural networks).
# Does so by dispatching into the `plot_fitted_funcs` function.
function plot_fitted_func!(rs::ReactionSystem, opt_sol, md::MeasuredData; kwargs...)
    plot_fitted_funcs!(rs, (u = [opt_sol],), md; kwargs...)
end

# Plots the fitted functions (potentially neural networks), and compares these with the true
# function. Does not handle functions of more than one variable. `f_idx` is the index of the
# function as stored in the `ReactionSystem` metadata. For systems with a single function (most)
# this will always be 1.
function plot_fitted_funcs!(rs::ReactionSystem, opt_sols, md::MeasuredData; f_idx = 1,
        range_expansion = 1.2, n_plt_vals = 500, plot_data_dist = false, include = :all,
        color_good_bounds = :green, color_good_outside = :darkolivegreen1, color_ok_bounds = :red1,
        color_ok_outside = :salmon1, color_bad_bounds = :grey88, color_bad_outside = :gray60, plot_outbounds = false,
        lw_bounds = 5, lw_outside = 2, la_bounds = 0.8, max_num = 200, real_x_range = true, linestyle_bounds = :solid,
        plot_true_rate = true, func_label = "Function rate", lw_true = 11, la_true = 0.5,
        plt_range = nothing, kwargs...)
    # Finds the function's input species.
    sp_in = only(get_dependent_sps(rs)[f_idx])
    (opt_sols isa Vector{ExpandedOptimizationSolution}) || (opt_sols = opt_sols.u)
    (include == :ok) && (opt_sols = filter(os -> os.objective < (1.5 ^ (true_loss(md) > 1 ? 1 : -1))*true_loss(md), opt_sols))
    (include == :good) && (opt_sols = filter(os -> os.objective < true_loss(md), opt_sols))

    # Error check.
    (length(get_dependent_sps(rs)[f_idx]) != 1)  && error("Currently, `plot_fitted_funcs` only supports functions of a single variable.")

    # Find the plotting range. Typically, we limit ourselves to values naturally occurring ground truth simulations.
    f_input_vals = if real_x_range
        [
            max.(sol[sp_in], 0.0) for sol in md.sols_true
            for opt_sol in opt_sols # Create one identical vector for each solution we wish to plot.
        ]
    else
        [
            vcat([max.(get_fitted_osol(opt_sol; cond)[sp_in],0.0) for cond in eachindex(md.u0s)]...)
            for opt_sol in opt_sols
        ]
    end
    f_min_ins = minimum.(f_input_vals)
    f_max_ins = maximum.(f_input_vals)
    plt_min, plt_max = if isnothing(plt_range)
        (minimum(f_min_ins) / range_expansion, maximum(f_max_ins) * range_expansion)
    else
        plt_range
    end
    plt_vals = range(plt_min; stop = plt_max, length = n_plt_vals)

    # Computes the expression corresponding to the true function (potentially inserts parameter values).
    true_rate = get_rates_of_interest(md.rs)[f_idx]
    subs_dict = Catalyst.symmap_to_varmap(md.rs, convert(Dict{Symbol,Float64}, relevant_ps(md.rs, md.ps_true)))
    true_rate = MTK.substitute(true_rate, subs_dict)

    # Compute the fitted function's value for all fits (different approach for nns and non-nns).
    fs_evaluated = if f_idx_is_nn(rs, f_idx)
        oprob_base, _, _, _, set_ps, _, _, _, _, _, _ = get_opt_prob_params(md, rs)
        [eval_nn(plt_vals, rs, f_idx_to_nn_idx(rs, f_idx), opt_sol.u, md, oprob_base, set_ps) for opt_sol in opt_sols]
    else
        symbolic_rate = get_rates_of_interest(rs)[f_idx]
        num_vals = [MTK.wrap(MTK.substitute(symbolic_rate, Dict(Pair.(fitted_parameters(rs), best_u(opt_sol)))))
            for opt_sol in opt_sols]
        [
            [
                MTK.value(MTK.substitute(num_val, Dict(getproperty(rs, sp_in) => sp_val)))
                for sp_val in plt_vals
            ]
            for num_val in num_vals
        ]
    end

    # Plots the fitted functions.
    for (idx, f_evaluated) in enumerate(fs_evaluated[1:min(max_num, length(fs_evaluated))])
        label = (idx == 1) ? func_label : ""
        good_fit = opt_sols[idx].objective < true_loss(md)
        ok_fit = opt_sols[idx].objective < (1.5 ^ (true_loss(md) > 1 ? 1 : -1))*true_loss(md)
        color1 = good_fit ? color_good_outside : (ok_fit ? color_ok_outside : color_bad_outside)
        color2 = good_fit ? color_good_bounds : (ok_fit ? color_ok_bounds : color_bad_bounds)
        la = [b ? la_bounds : plot_outbounds ? 0.1 : 0.0 for b in (f_min_ins[idx] .<= plt_vals .<= f_max_ins[idx])]
        lw = [b ? lw_bounds : lw_outside for b in (f_min_ins[idx] .<= plt_vals .<= f_max_ins[idx])]
        color = [b ? color2 :  color1 for b in (f_min_ins[idx] .<= plt_vals .<= f_max_ins[idx])]
        linestyle = [b ? linestyle_bounds : :dot for b in (f_min_ins[idx] .<= plt_vals .<= f_max_ins[idx])]
        if !plot_outbounds
            la = la_bounds
            lw = lw_bounds
            color = color2
            linestyle = linestyle_bounds
        end
        plot!(plt_vals, f_evaluated; label, lw, la, color, linestyle)
    end

    # Plots the true rate's value across the target range.
    true_vals = [MTK.substitute(Symbolics.unwrap(true_rate), Catalyst.symmap_to_varmap(rs, [sp_in => val])) for val in plt_vals]
    plot_true_rate && plot!(plt_vals, true_vals; label = "True rate", lw = lw_true, la = la_true, color = :darkorange, linestyle = :dash)
    plt = plot!()
    ylimit = (max(ylims(plt)[1], 0.0), min(ylims(plt)[2], 1.3*maximum(true_vals)))

    # Plots information of the range of values that is covered by the data. This includes potentially
    # plotting the distribution of data values across the range.
    if plot_data_dist
        plot!([data_range[1], data_range[1]], collect(ylims(plt)); ylimit, color = :black, la = 0.4,
            lw = 2, linestyle = :dot, label = "",)
        plot!([data_range[2], data_range[2]], collect(ylims(plt)); ylimit, color = :black, la = 0.4,
            lw = 2, linestyle = :dot, label = "")
        val_density = kde(all_vals)
        rescaled_density = ylimit[2] * val_density.density / maximum(val_density.density)
        plot!(val_density.x, rescaled_density; ylimit, color = :black, lw = 2, la = 0.6, linestyle = :dot, label = "")
    end

    # Makes final adjustments and returns the plot.
    plot!(; xlims = (plt_min, plt_max), xguide = sp_in, kwargs...)
end

# Dispatch for when a `SavedRun` object is provided.
function plot_fitted_funcs!(sr::SavedRun, model_idx::Int64; f_idx = 1, range_expansion = 1.2, n_plt_vals = 500,
        plot_data_dist = false, include = :all, color = model_type_color(sr.models[model_idx]), kwargs...)
    rs = sr.models[model_idx]
    opt_sols = sr.opt_sols[model_idx]
    md = sr.md
    return plot_fitted_funcs!(rs, opt_sols, md; f_idx, color_ok_bounds = color, color_good_bounds = color, range_expansion, n_plt_vals, plot_data_dist, include, kwargs...)
end

# Plots the loss function evaluation traces across all optimisation runs.
function plot_loss_traces(opt_sol, md; ymax = 1e5)
    plot(xlimit = (1, maximum(length(sol.trace) for sol in opt_sol)))
    for sol in opt_sol
        color = sol.objective < true_loss(md) ? :gold :
            (sol.objective < 1.5*true_loss(md) ? :darkorange2 : :grey)
        plot!(sol.min_trace[1:end]; label = "", lw = 3, la = 0.7, color)
        scatter!((sol.opt_sol1.stats.iterations, sol.trace[sol.opt_sol1.stats.iterations]);
            markersize = 3, markerstrokewidth = 1, label = "", color)
    end
    plot!(yaxis = :log10, ylimit = (-Inf, ymax))
end

### Defines Plotting Constants ###

# Function for retrieving a colour for a specific index.
get_color(idx) = (idx <= length(plotting_colors)) ? plotting_colors[idx] : [idx - length(plotting_colors)]
get_color_alt(idx) = (idx <= length(plotting_colors_alt)) ? plotting_colors_alt[idx] : [idx - length(plotting_colors_alt)]

# The colours for the plots. Each triplet correspond to colors for: (1) True simulation,
# (2) Measured values, and (3) Fitted simulation.
plotting_colors = (
    (:lightskyblue1, :deepskyblue, :blue1),
    (:coral, :coral2, :crimson),
    (:lawngreen, :chartreuse2, :green4),
    (:mediumpurple1, :orchid1, :darkviolet),
    (:lightgoldenrod1, :gold, :darkgoldenrod1),
    (:cyan1, :turquoise1, :turquoise3)
)

# An alternative scheme, which uses different colors for the 3 models.
plotting_colors_alt = (
    (:lightskyblue1, :lightskyblue1, :dodgerblue3, :dodgerblue, :slateblue2),
    (:coral, :coral, :red2, :firebrick1, :crimson),
    (:lawngreen, :lawngreen, :springgreen2, :chartreuse3, :olivedrab),
    (:mediumpurple1, :mediumpurple1, :orchid1, :darkorchid2, :purple),
    (:lightgoldenrod1, :lightgoldenrod1, :goldenrod1, :khaki3, :darkgoldenrod3),
    (:cyan1, :cyan1, :turquoise1, :turquoise3, :cyan4)
)

# An alternative scheme, which uses different colors for the 3 models.
plotting_colors_ic = (
    (:lightskyblue1, :lightskyblue1, :dodgerblue, :cadetblue, :deepskyblue3, :slateblue2),
    (:coral, :coral,  :firebrick1, :tomato, :indianred,:crimson),
    (:lawngreen, :lawngreen, :chartreuse3, :lightseagreen, :seagreen, :olivedrab),
)

# An alternative scheme, which uses different colors for the 4 models.
plotting_colors_linpathway = (
    (:lightskyblue1, :lightskyblue1, :dodgerblue, :cadetblue, :deepskyblue3, :slateblue2),
    (:coral, :coral,  :firebrick1, :tomato, :indianred,:crimson),
    (:lawngreen, :lawngreen, :chartreuse3, :lightseagreen, :seagreen, :olivedrab),
    (:mediumpurple1, :mediumpurple1, :orchid1, :darkorchid2, :purple3, :mediumpurple),
)
