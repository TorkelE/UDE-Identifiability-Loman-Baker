### Preparations ###

# Sets the correct directory. Fetches packages.
begin
    cd("$(@__DIR__)/../..")
    using Pkg; Pkg.activate(".")
    using NeuralCRNsIdentifiability
    using OrdinaryDiffEq
    using Plots
    using Serialization
    import Plots: mm
end

# Set plotting defaults.
default(framestyle = :box)

# Fetches the investigations.
@time begin
    # Fetch raw data fields (full scans takes a lot of space).
    mosquitoes_data = deserialize("Scans/mosquitoes_data_variation_scan/extracted_vals.jls")
    extended_sal_data = deserialize("Scans/extended_sal_data_variation_scan/extracted_vals.jls")
end

### Figure generation ###

# Function for plotting an extracted scan.
function plot_scan(scan_evaluation, property_idx, model_idx; xguide = "Data quality", kwargs...)
    means = reverse(scan_evaluation[1][model_idx][property_idx])
    stds = reverse(scan_evaluation[2][model_idx][property_idx])
    plot(eachindex(means), means; ribbon = stds, xguide, xlimit = (1, length(means)), kwargs...)
end

# Plots the full grid of all evaluations for the mosquito model.
begin
    # Fetches the correct colors.
    c1 = :darkorange
    c2 = :purple
    c3 = :green

    # Creates the likelihood profile plots.
    lw = 7
    lp_plt1 = plot_scan(mosquitoes_data, 1, 1; ylimit = (0.0, 0.02), yticks = [0.0, 0.005, 0.01, 0.015, 0.02], color = c1, lw)
    lp_plt2 = plot_scan(mosquitoes_data, 1, 2; ylimit = (0.0, 0.2), yticks = [0.0, 0.05, 0.1, 0.15, 0.2], color = c2, lw)
    lp_plt3 = plot_scan(mosquitoes_data, 1, 3; ylimit = (0.0, 6.0), yticks = [0.0, 1.5, 3.0, 4.5, 6.0], color = c3, lw)

    # Creates the functional mean diff plots.
    lw = 7
    funcdiff_plt1 = plot_scan(mosquitoes_data, 3, 1; ylimit = (-0.05, 1.05), yticks = [0.0, 0.25, 0.5, 0.75, 1.0], color = c1, lw)
    funcdiff_plt2 = plot_scan(mosquitoes_data, 3, 2; ylimit = (0.0, 4.0), yticks = [0, 1, 2, 3, 4], color = c2, lw)
    funcdiff_plt3 = plot_scan(mosquitoes_data, 3, 3; ylimit = (0.0, 60.0), yticks = [0, 15, 30, 45, 60], color = c3, lw)

    # Creates the prediction error plots.
    lw = 7
    prederror_plt1 = plot_scan(mosquitoes_data, 2, 1; ylimit = (0.0, 2.4), yticks = [0.0, 0.8, 1.6, 2.4], color = c1, lw)
    prederror_plt2 = plot_scan(mosquitoes_data, 2, 2; ylimit = (0.0, 15.0), yticks = [0, 5, 10, 15], color = c2, lw)
    prederror_plt3 = plot_scan(mosquitoes_data, 2, 3; ylimit = (0.0, 36.0), yticks = [0, 12, 24, 36], color = c3, lw)

    # Makes the final combined plot
    xticks = [1, 4, 7, 10]
    xguide = ""
    bottom_margin = [-2mm -2mm -2mm -2mm -2mm -2mm -2mm -2mm -2mm ]
    left_margin = [-1mm -1mm -1mm -1mm -1mm -1mm -1mm -1mm -1mm]
    mosquitoes_plt = plot(
        lp_plt1, lp_plt2, lp_plt3, funcdiff_plt1, funcdiff_plt2, funcdiff_plt3, prederror_plt1, prederror_plt2, prederror_plt3;
        xguide, xticks, bottom_margin, left_margin, layout = (3,3), size = (1350, 850)
    )
    make_fig(mosquitoes_plt, "Mosquitoes/FigureComponents/mosquitoes_datavar_base")
end

# Plots the full grid of all evaluations for the extended self-activation loop model.
begin
    # Fetches the correct colors.
    c1 = :darkorange
    c2 = :purple
    c3 = :green

    # Creates the likelihood profile plots.
    lw = 7
    lp_plt1 = plot_scan(extended_sal_data, 1, 1; ylimit = (0.0, 0.01), yticks = [0.0, 0.0025, 0.005, 0.0075, 0.01], color = c1, lw)
    lp_plt2 = plot_scan(extended_sal_data, 1, 2; ylimit = (0.0, 1.0), yticks = [0.0, 0.25, 0.5, 0.75, 1.0], color = c2, lw)
    lp_plt3 = plot_scan(extended_sal_data, 1, 3; ylimit = (0.0, 1.0), yticks = [0.0, 0.25, 0.5, 0.75, 1.0], color = c3, lw)

    # Creates the functional mean diff plots.
    lw = 7
    funcdiff_plt1 = plot_scan(extended_sal_data, 3, 1; ylimit = (-0.05, 1.05), yticks = [0.0, 0.25, 0.5, 0.75, 1.0], color = c1, lw)
    funcdiff_plt2 = plot_scan(extended_sal_data, 3, 2; ylimit = (0.0, 60.0), yticks = [0, 15, 30, 45, 60], color = c2, lw)
    funcdiff_plt3 = plot_scan(extended_sal_data, 3, 3; ylimit = (0.0, 60.0), yticks = [0, 15, 30, 45, 60], color = c3, lw)

    # Creates the prediction error plots.
    lw = 7
    prederror_plt1 = plot_scan(extended_sal_data, 2, 1; ylimit = (0.0, 3.6), yticks = [0, 1.2, 2.4, 3.6], color = c1, lw)
    prederror_plt2 = plot_scan(extended_sal_data, 2, 2; ylimit = (0.0, 36.0), yticks = [0, 12, 24, 36], color = c2, lw)
    prederror_plt3 = plot_scan(extended_sal_data, 2, 3; ylimit = (0.0, 36.0), yticks = [0, 12, 24, 36], color = c3, lw)

    # Makes the final combined plot
    xticks = [1, 4, 7, 10]
    xguide = ""
    bottom_margin = [-2mm -2mm -2mm -2mm -2mm -2mm -2mm -2mm -2mm ]
    left_margin = [-1mm -1mm -1mm -1mm -1mm -1mm -1mm -1mm -1mm]
    extended_sal_plt = plot(
        lp_plt1, lp_plt2, lp_plt3, funcdiff_plt1, funcdiff_plt2, funcdiff_plt3, prederror_plt1, prederror_plt2, prederror_plt3;
        xguide, xticks, bottom_margin, left_margin, layout = (3,3), size = (1350, 850)
    )
    make_fig(extended_sal_plt, "Mosquitoes/FigureComponents/extended_sal_datavar_base")
end

# Mousqitoes sample points.
begin
    # Plots examples cases (mosquitoes).
    run_tag = "mosquitoes_data_samples_scan_1"
    run1 = sampled_fits1.runs[1]
    run2 = sampled_fits2.runs[8]
    run3 = sampled_fits3.runs[17]
    run4 = sampled_fits4.runs[23]

    # Plots the fits.
    ylimit = (0.0, 1.5); 1; 1; xguide = ""; lws = [6, 5, 4]; ms_true = 5; lw_true = 8; la_true = 1.0; las = [0.7, 0.6, 0.5]; legend = :none;
    sim_plt_1 = plot_fitted_sims_combined(run1; ylimit, yticks = [0.0, 0.5, 1.0, 1.5], xticks = [0.0, 20.0, 40.0, 60.0, 80.0], xguide, lws, ms_true, lw_true, la_true, las)
    sim_plt_2 = plot_fitted_sims_combined(run2; ylimit, yticks = [0.0, 0.5, 1.0, 1.5], xticks = [0.0, 20.0, 40.0, 60.0, 80.0], xguide, lws, ms_true, lw_true, la_true, las, legend)
    sim_plt_3 = plot_fitted_sims_combined(run3; ylimit, yticks = [0.0, 0.5, 1.0, 1.5], xticks = [0.0, 20.0, 40.0, 60.0, 80.0], xguide, lws, ms_true, lw_true, la_true, las, legend)
    sim_plt_4 = plot_fitted_sims_combined(run4; ylimit, yticks = [0.0, 0.5, 1.0, 1.5], xticks = [0.0, 20.0, 40.0, 60.0, 80.0], xguide, lws, ms_true, lw_true, la_true, las, legend)

    # Plots the profiles.
    lws = [8, 6, 5]; las = [0.8, 0.7, 0.6]; no_marker = true; ymin = -4.2; xguide = ""; yguide = ""; xlimit = (-2.0, 4.0);
    l_plt_1 = plot_combined_lprofiles(run1, :r; lws, xlimit, no_marker, ymin, xguide, yguide, yticks = [-4, -2, 0], xticks = [-2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0], las)
    l_plt_2 = plot_combined_lprofiles(run2, :r; lws, xlimit, no_marker, ymin, xguide, yguide, yticks = [-4, -2, 0], xticks = [-2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0], las)
    l_plt_3 = plot_combined_lprofiles(run3, :r; lws, xlimit, no_marker, ymin, xguide, yguide, yticks = [-4, -2, 0], xticks = [-2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0], las)
    l_plt_4 = plot_combined_lprofiles(run4, :r; lws, xlimit, no_marker, ymin, xguide, yguide, yticks = [-4, -2, 0], xticks = [-2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0], las)

    # Plots the fitted functions.
    xlimit = (0.0, 0.75); ylimit = (0.0, 1.0); xguide = ""; lw_true = 13; la_true = 0.6; legend = :none;
    plt_funcs_1 = plot_combined_fitted_funcs(run1; xlimit, ylimit, xticks = [0.0, 0.25, 0.5, 0.75], yticks = [0.0, 0.5, 1.0], xguide, lw_true, la_true, legend)
    plt_funcs_2 = plot_combined_fitted_funcs(run2; xlimit, ylimit, xticks = [0.0, 0.25, 0.5, 0.75], yticks = [0.0, 0.5, 1.0], xguide, lw_true, la_true, legend)
    plt_funcs_3 = plot_combined_fitted_funcs(run3; xlimit, ylimit = (0.0, 2.0), xticks = [0.0, 0.25, 0.5, 0.75], yticks = [0.0, 1.0, 2.0], xguide, lw_true, la_true, legend)
    plt_funcs_4 = plot_combined_fitted_funcs(run4; xlimit, ylimit = (0.0, 4.0), xticks = [0.0, 0.25, 0.5, 0.75], yticks = [0.0, 1.0, 2.0, 3.0, 4.0], xguide, lw_true, la_true, legend)

    # Makes the final combined plot
    bottom_margin = [28mm 28mm 28mm 28mm 20mm 20mm 20mm 20mm 0mm 0mm 0mm 0mm]
    mosquitoes_samples_plt = plot(sim_plt_4, sim_plt_3, sim_plt_2, sim_plt_1, l_plt_4, l_plt_3, l_plt_2, l_plt_1, plt_funcs_4, plt_funcs_3, plt_funcs_2, plt_funcs_1;
        layout = (3,4), size = (1800, 1200), legendfontsize = 3, left_margin = -5mm, bottom_margin)
    make_fig(mosquitoes_samples_plt, "Mosquitoes/FigureComponents/mosquitoes_sample_fits_base")
end

# Extended activation loop sample points.
begin
    # Plots examples cases (extended SAL).
    selected_runs = vcat(1:40, 3*40 .+ (1:40), 6*40 .+ (1:40), 9*40 .+ (1:40))
    selected_runs[1] = selected_runs[2]
    selected_runs[143:end] .= selected_runs[142]
    run_tag = "extended_sal_data_variation_scan"
    sampled_fits = NeuralCRNsIdentifiability.SavedScan(run_tag, "Scans/$run_tag"; selected_runs)
    idxs = [24, 80, 109, 123] # Select ones where profile likelihoods are not that jagged (computing these for UDEs in a pretty way is hard).

    # Plots the fits.
    ylimit = (0.0, 2.1); yticks = [0.0, 0.7, 1.4, 2.1]; xticks = [0.0, 5.0, 10.0, 15.0, 20.0]; xguide = ""; lws = [6, 5, 4]; ms_true = 5; lw_true = 8; la_true = 1.0; las = [0.7, 0.6, 0.5]; legend = :none;
    sim_plt_1 = plot_fitted_sims_combined(sampled_fits.runs[idxs[1]]; ylimit, yticks, xticks, xguide, lws, ms_true, lw_true, la_true, las)
    sim_plt_2 = plot_fitted_sims_combined(sampled_fits.runs[idxs[2]]; ylimit, yticks, xticks, xguide, lws, ms_true, lw_true, la_true, las, legend)
    sim_plt_3 = plot_fitted_sims_combined(sampled_fits.runs[idxs[3]]; ylimit, yticks, xticks, xguide, lws, ms_true, lw_true, la_true, las, legend)
    sim_plt_4 = plot_fitted_sims_combined(sampled_fits.runs[idxs[4]]; ylimit, yticks, xticks, xguide, lws, ms_true, lw_true, la_true, las, legend)

    # Plots the profiles.
    lws = [8, 6, 5]; las = [0.8, 0.7, 0.6]; no_marker = true; ymin = -4.2; xguide = ""; yguide = ""; yticks = [-4, -2, 0]; xlimit = (-2.0, 1.0); xticks = [-2.0, -1.0, 0.0, 1.0];
    l_plt_1 = plot_combined_lprofiles(sampled_fits.runs[idxs[1]], :d; lws, xlimit, no_marker, ymin, xguide, yguide, yticks, las)
    l_plt_2 = plot_combined_lprofiles(sampled_fits.runs[idxs[2]], :d; lws, xlimit, no_marker, ymin, xguide, yguide, yticks, las)
    l_plt_3 = plot_combined_lprofiles(sampled_fits.runs[idxs[3]], :d; lws, xlimit, no_marker, ymin, xguide, yguide, yticks, las)
    l_plt_4 = plot_combined_lprofiles(sampled_fits.runs[idxs[4]], :d; lws, xlimit, no_marker, ymin, xguide, yguide, yticks, las)

    # Plots the fitted functions.
    xlimit = (0.0, 2.0); ylimit = (0.0, 1.5); xticks = [0.0, 1.0, 2.0]; yticks = [0.0, 0.75, 1.5]; xguide = ""; lw_true = 13; la_true = 0.6; legend = :none;
    plt_funcs_1 = plot_combined_fitted_funcs(sampled_fits.runs[idxs[1]]; xlimit, ylimit, xticks, yticks, xguide, lw_true, la_true, legend)
    plt_funcs_2 = plot_combined_fitted_funcs(sampled_fits.runs[idxs[2]]; xlimit, ylimit, xticks, yticks, xguide, lw_true, la_true, legend)
    plt_funcs_3 = plot_combined_fitted_funcs(sampled_fits.runs[idxs[3]]; xlimit, ylimit, xticks, yticks, xguide, lw_true, la_true, legend)
    plt_funcs_4 = plot_combined_fitted_funcs(sampled_fits.runs[idxs[4]]; xlimit, ylimit = (0.0, 4.0), xticks, yticks = [0.0, 1.0, 2.0, 3.0, 4.0], xguide, lw_true, la_true, legend)

    # Makes the final combined plot
    bottom_margin = [28mm 28mm 28mm 28mm 20mm 20mm 20mm 20mm 0mm 0mm 0mm 0mm]
    extended_sal_samples_plt = plot(sim_plt_4, sim_plt_3, sim_plt_2, sim_plt_1, l_plt_4, l_plt_3, l_plt_2, l_plt_1, plt_funcs_4, plt_funcs_3, plt_funcs_2, plt_funcs_1;
        layout = (3,4), size = (1800, 1200), legendfontsize = 3, left_margin = -5mm, bottom_margin)
    make_fig(extended_sal_samples_plt, "Mosquitoes/FigureComponents/extended_sal_sample_fits_base")
end

### Check Plots ###

sampled_fits1 = NeuralCRNsIdentifiability.SavedScan(run_tag, "Scans/mosquitoes_data_samples_scan_1")
sampled_fits2 = NeuralCRNsIdentifiability.SavedScan(run_tag, "Scans/mosquitoes_data_samples_scan_2")
sampled_fits3 = NeuralCRNsIdentifiability.SavedScan(run_tag, "Scans/mosquitoes_data_samples_scan_3")
sampled_fits4 = NeuralCRNsIdentifiability.SavedScan(run_tag, "Scans/mosquitoes_data_samples_scan_4")
