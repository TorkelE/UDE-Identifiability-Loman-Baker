### Preparations ###

# Sets the correct directory. Fetches packages.
begin
    cd("$(@__DIR__)/../..")
    using Pkg; Pkg.activate(".")
    using NeuralCRNsIdentifiability
    using Plots
    import Plots: mm
end

# Set plotting defaults.
default(framestyle = :box)

# Fetches the investigations.
begin
    # Fetches the main investigation (repeats included in here).
    run_tag_main = "extended_sal_basic_scan"
    check_scan(run_tag_main; print_info = true)
    investigation_main = load_investigation(run_tag_main)

    # Fetches the first alternative parameters investigation.
    run_tag_alt_ps = "extended_sal_basic_scan_alt_parameters"
    check_scan(run_tag_alt_ps; print_info = true)
    investigation_alt_ps = load_investigation(run_tag_alt_ps)

    # Fetches the second (low concentration) alternative parameters investigation.
    run_tag_alt_ps_low_conc = "extended_sal_basic_scan_alt_parameters_low_conc"
    check_scan(run_tag_alt_ps_low_conc; print_info = true)
    investigation_alt_ps_low_conc = load_investigation(run_tag_alt_ps_low_conc)

    # Fetches the alternative model (GWO) investigation.
    run_tag_gwo = "gwo_basic_scan"
    check_scan(run_tag_gwo; print_info = true)
    investigation_gwo = load_investigation(run_tag_gwo)

    # Fetches the extensive data fitting run.
    run_tag_high_data = "extended_sal_extensive_data_run"
    investigation_high_data = load_investigation(run_tag_high_data)
end


### Figure generation - Main ###

# Investigation of practical identifiability for the extended self-activation loop.
begin
    # Generates the simulation plots (for each data set).
    ylimit = (0.0, 6.0); yticks = [0.0, 2.0, 4.0, 6.0]; xticks = [0.0, 15.0, 30.0, 45.0]; xguide = ""; lws = [10, 8, 6]; ms_true = 7; lw_true = 12; la_true = 1.0; las = [0.8, 0.7, 0.6]; legend = :none;
    sim_plt_1 = plot_fitted_sims_combined(investigation_main.runs[1]; ylimit, yticks, xticks, xguide, lws, ms_true, lw_true, la_true, las, legend)
    sim_plt_2 = plot_fitted_sims_combined(investigation_main.runs[2]; ylimit, yticks, xticks, xguide, lws, ms_true, lw_true, la_true, las, legend)
    sim_plt_3 = plot_fitted_sims_combined(investigation_main.runs[3]; ylimit, yticks, xticks, xguide, lws, ms_true, lw_true, la_true, las, legend)

    # Generates the likelihood profile plots (for each data set).
    no_marker = true; ymin = -4.2; xguide = ""; yguide = ""; yticks = [-4, -2, 0]; xlimit = (-1.0, 1.0); xticks = [-1.0, 0.0, 1.0];
    pl_plt_1 = plot_combined_lprofiles(investigation_main.runs[1], :d; lws = [10, 10, 10], no_marker, ymin, xguide, yguide, yticks, las = [0.8, 0.8, 0.8], xlimit, xticks)
    pl_plt_2 = plot_combined_lprofiles(investigation_main.runs[2], :d; lws = [10, 8, 5], no_marker, ymin, xguide, yguide, yticks, las = [0.8, 0.7, 0.6], xlimit, xticks)
    pl_plt_3 = plot_combined_lprofiles(investigation_main.runs[3], :d; lws = [10, 8, 5], no_marker, ymin, xguide, yguide, yticks, las = [0.8, 0.7, 0.6], xlimit, xticks)

    # Generates the ensemble function plots (for each data set).
    xlimit = (0.0, 4.5); ylimit = (0.0, 1.2); xticks = [0.0, 1.5, 3.0, 4.5]; yticks = [0.0, 0.4, 0.8, 1.2]; xguide = ""; lw_true = 13; la_true = 0.6;
    plt_funcs_1 = plot_combined_fitted_funcs(investigation_main.runs[1]; plt_range = (0.0, 4.5), ylimit, xticks, xlimit, yticks, xguide, lw_true, la_true, max_num = 200, max_num_nn = 75)
    plt_funcs_2 = plot_combined_fitted_funcs(investigation_main.runs[2]; plt_range = (0.0, 4.5), ylimit, xticks, xlimit, yticks, xguide, lw_true, la_true, max_num = 30)
    plt_funcs_3 = plot_combined_fitted_funcs(investigation_main.runs[3]; plt_range = (0.0, 4.5), ylimit, xticks, xlimit, yticks, xguide, lw_true, la_true, max_num = 30)

    # Generates the combined plot.
    bottom_margin = [28mm 28mm 28mm 20mm 20mm 20mm 0mm 0mm 0mm]
    extended_sal_main_plt = plot(sim_plt_1, sim_plt_2, sim_plt_3, pl_plt_1, pl_plt_2, pl_plt_3, plt_funcs_1, plt_funcs_2, plt_funcs_3;
        layout = (3,3), size = (1350, 1200), legendfontsize = 3, left_margin = -5mm, bottom_margin)
    make_fig(extended_sal_main_plt, "Extended_SA_Loop/FigureComponents/extended_sal_pract_ident_base")
end

### Figure generation - Supplementary ###

# Main figure - repeat for additional data draw (1).
let
    # Generates the simulation plots (for each data set).
    ylimit = (0.0, 6.0); yticks = [0.0, 2.0, 4.0, 6.0]; xticks = [0.0, 15.0, 30.0, 45.0]; xguide = ""; lws = [10, 8, 6]; ms_true = 7; lw_true = 12; la_true = 1.0; las = [0.8, 0.7, 0.6]; legend = :none;
    sim_plt_1 = plot_fitted_sims_combined(investigation_main.runs[4]; ylimit, yticks, xticks, xguide, lws, ms_true, lw_true, la_true, las, legend)
    sim_plt_2 = plot_fitted_sims_combined(investigation_main.runs[5]; ylimit, yticks, xticks, xguide, lws, ms_true, lw_true, la_true, las, legend)
    sim_plt_3 = plot_fitted_sims_combined(investigation_main.runs[6]; ylimit, yticks, xticks, xguide, lws, ms_true, lw_true, la_true, las, legend)

    # Generates the likelihood profile plots (for each data set).
    no_marker = true; ymin = -4.2; xguide = ""; yguide = ""; yticks = [-4, -2, 0]; xlimit = (-1.0, 1.0); xticks = [-1.0, 0.0, 1.0];
    pl_plt_1 = plot_combined_lprofiles(investigation_main.runs[4], :d; lws = [10, 10, 10], no_marker, ymin, xguide, yguide, yticks, las = [0.8, 0.8, 0.8], xlimit, xticks)
    pl_plt_2 = plot_combined_lprofiles(investigation_main.runs[5], :d; lws = [10, 8, 5], no_marker, ymin, xguide, yguide, yticks, las = [0.8, 0.7, 0.6], xlimit, xticks)
    pl_plt_3 = plot_combined_lprofiles(investigation_main.runs[6], :d; lws = [10, 8, 5], no_marker, ymin, xguide, yguide, yticks, las = [0.8, 0.7, 0.6], xlimit, xticks)

    # Generates the ensemble function plots (for each data set).
    xlimit = (0.0, 4.5); ylimit = (0.0, 1.2); xticks = [0.0, 1.5, 3.0, 4.5]; yticks = [0.0, 0.4, 0.8, 1.2]; xguide = ""; lw_true = 13; la_true = 0.6;
    plt_funcs_1 = plot_combined_fitted_funcs(investigation_main.runs[4]; plt_range = (0.0, 4.5), ylimit, xticks, xlimit, yticks, xguide, lw_true, la_true)
    plt_funcs_2 = plot_combined_fitted_funcs(investigation_main.runs[5]; plt_range = (0.0, 4.5), ylimit, xticks, xlimit, yticks, xguide, lw_true, la_true, max_num = 30)
    plt_funcs_3 = plot_combined_fitted_funcs(investigation_main.runs[6]; plt_range = (0.0, 4.5), ylimit, xticks, xlimit, yticks, xguide, lw_true, la_true, max_num = 30)

    # Generates the combined plot.
    bottom_margin = [28mm 28mm 28mm 20mm 20mm 20mm 0mm 0mm 0mm]
    extended_sal_rep1_plt = plot(sim_plt_1, sim_plt_2, sim_plt_3, pl_plt_1, pl_plt_2, pl_plt_3, plt_funcs_1, plt_funcs_2, plt_funcs_3;
        layout = (3,3), size = (1350, 1200), legendfontsize = 3, left_margin = -5mm, bottom_margin)
    make_fig(extended_sal_rep1_plt, "Extended_SA_Loop/FigureComponents/extended_sal_pract_ident_rep1_base")
end

# Main figure - repeat for additional data draw (2).
let
    # Generates the simulation plots (for each data set).
    ylimit = (0.0, 6.0); yticks = [0.0, 2.0, 4.0, 6.0]; xticks = [0.0, 15.0, 30.0, 45.0]; xguide = ""; lws = [10, 8, 6]; ms_true = 7; lw_true = 12; la_true = 1.0; las = [0.8, 0.7, 0.6]; legend = :none;
    sim_plt_1 = plot_fitted_sims_combined(investigation_main.runs[7]; ylimit, yticks, xticks, xguide, lws, ms_true, lw_true, la_true, las, legend)
    sim_plt_2 = plot_fitted_sims_combined(investigation_main.runs[8]; ylimit, yticks, xticks, xguide, lws, ms_true, lw_true, la_true, las, legend)
    sim_plt_3 = plot_fitted_sims_combined(investigation_main.runs[9]; ylimit, yticks, xticks, xguide, lws, ms_true, lw_true, la_true, las, legend)

    # Generates the likelihood profile plots (for each data set).
    no_marker = true; ymin = -4.2; xguide = ""; yguide = ""; yticks = [-4, -2, 0]; xlimit = (-1.0, 1.0); xticks = [-1.0, 0.0, 1.0];
    pl_plt_1 = plot_combined_lprofiles(investigation_main.runs[7], :d; lws = [10, 10, 10], no_marker, ymin, xguide, yguide, yticks, las = [0.8, 0.8, 0.8], xlimit, xticks)
    pl_plt_2 = plot_combined_lprofiles(investigation_main.runs[8], :d; lws = [10, 8, 5], no_marker, ymin, xguide, yguide, yticks, las = [0.8, 0.7, 0.6], xlimit, xticks)
    pl_plt_3 = plot_combined_lprofiles(investigation_main.runs[9], :d; lws = [10, 8, 5], no_marker, ymin, xguide, yguide, yticks, las = [0.8, 0.7, 0.6], xlimit, xticks)

    # Generates the ensemble function plots (for each data set).
    xlimit = (0.0, 4.5); ylimit = (0.0, 1.2); xticks = [0.0, 1.5, 3.0, 4.5]; yticks = [0.0, 0.4, 0.8, 1.2]; xguide = ""; lw_true = 13; la_true = 0.6;
    plt_funcs_1 = plot_combined_fitted_funcs(investigation_main.runs[7]; plt_range = (0.0, 4.5), ylimit, xticks, xlimit, yticks, xguide, lw_true, la_true)
    plt_funcs_2 = plot_combined_fitted_funcs(investigation_main.runs[8]; plt_range = (0.0, 4.5), ylimit, xticks, xlimit, yticks, xguide, lw_true, la_true, max_num = 30)
    plt_funcs_3 = plot_combined_fitted_funcs(investigation_main.runs[9]; plt_range = (0.0, 4.5), ylimit, xticks, xlimit, yticks, xguide, lw_true, la_true, max_num = 30)

    # Generates the combined plot.
    bottom_margin = [28mm 28mm 28mm 20mm 20mm 20mm 0mm 0mm 0mm]
    extended_sal_rep2_plt = plot(sim_plt_1, sim_plt_2, sim_plt_3, pl_plt_1, pl_plt_2, pl_plt_3, plt_funcs_1, plt_funcs_2, plt_funcs_3;
        layout = (3,3), size = (1350, 1200), legendfontsize = 3, left_margin = -5mm, bottom_margin)
    make_fig(extended_sal_rep2_plt, "Extended_SA_Loop/FigureComponents/extended_sal_pract_ident_rep2_base")
end

# Main figure - repeat for an additional parameter set.
let
    # Generates the simulation plots (for each data set).
    ylimit = (0.0, 15.0); yticks = [0.0, 5.0, 10.0, 15.0]; xticks = [0.0, 15.0, 30.0, 45.0]; xguide = ""; lws = [10, 8, 6]; ms_true = 7; lw_true = 12; la_true = 1.0; las = [0.8, 0.7, 0.6]; legend = :none;
    sim_plt_1 = plot_fitted_sims_combined(investigation_alt_ps.runs[1]; ylimit, yticks, xticks, xguide, lws, ms_true, lw_true, la_true, las, legend)
    sim_plt_2 = plot_fitted_sims_combined(investigation_alt_ps.runs[2]; ylimit, yticks, xticks, xguide, lws, ms_true, lw_true, la_true, las, legend)
    sim_plt_3 = plot_fitted_sims_combined(investigation_alt_ps.runs[3]; ylimit, yticks, xticks, xguide, lws, ms_true, lw_true, la_true, las, legend)

    # Generates the likelihood profile plots (for each data set).
    no_marker = true; ymin = -4.2; xguide = ""; yguide = ""; yticks = [-4, -2, 0]; xlimit = (-1.0, 1.0); xticks = [-2.0, -1.0, 0.0, 1.0];
    pl_plt_1 = plot_combined_lprofiles(investigation_alt_ps.runs[1], :d; lws = [10, 10, 10], no_marker, ymin, xguide, yguide, yticks, las = [0.8, 0.8, 0.8], xlimit, xticks)
    pl_plt_2 = plot_combined_lprofiles(investigation_alt_ps.runs[2], :d; lws = [10, 8, 5], no_marker, ymin, xguide, yguide, yticks, las = [0.8, 0.7, 0.6], xlimit, xticks)
    pl_plt_3 = plot_combined_lprofiles(investigation_alt_ps.runs[3], :d; lws = [10, 8, 5], no_marker, ymin, xguide, yguide, yticks, las = [0.8, 0.7, 0.6], xlimit, xticks)

    # Generates the ensemble function plots (for each data set).
    xlimit = (0.0, 9.0); ylimit = (0.0, 2.7); xticks = [0.0, 3.0, 6.0, 9.0]; yticks = [0.0, 0.9, 1.8, 2.7]; xguide = ""; lw_true = 13; la_true = 0.6;
    plt_funcs_1 = plot_combined_fitted_funcs(investigation_alt_ps.runs[1]; plt_range = (0.0, 9.0), ylimit, xticks, xlimit, yticks, xguide, lw_true, la_true)
    plt_funcs_2 = plot_combined_fitted_funcs(investigation_alt_ps.runs[2]; plt_range = (0.0, 9.0), ylimit, xticks, xlimit, yticks, xguide, lw_true, la_true, max_num = 30)
    plt_funcs_3 = plot_combined_fitted_funcs(investigation_alt_ps.runs[3]; plt_range = (0.0, 9.0), ylimit, xticks, xlimit, yticks, xguide, lw_true, la_true, max_num = 30)

    # Generates the combined plot.
    bottom_margin = [28mm 28mm 28mm 20mm 20mm 20mm 0mm 0mm 0mm]
    extended_sal_alt_ps_plt = plot(sim_plt_1, sim_plt_2, sim_plt_3, pl_plt_1, pl_plt_2, pl_plt_3, plt_funcs_1, plt_funcs_2, plt_funcs_3;
        layout = (3,3), size = (1350, 1200), legendfontsize = 3, left_margin = -5mm, bottom_margin)
    make_fig(extended_sal_alt_ps_plt, "Extended_SA_Loop/FigureComponents/extended_sal_pract_ident_alt_ps_base")
end

# Main figure - repeat for an additional parameter set (but with distinct, low concentration, behaviours).
let
    # Generates the simulation plots (for each data set).
    ylimit = (0.0, 2.1); yticks = [0.0, 0.7, 1.4, 2.1]; xticks = [0.0, 15.0, 30.0, 45.0]; xguide = ""; lws = [10, 8, 6]; ms_true = 7; lw_true = 12; la_true = 1.0; las = [0.8, 0.7, 0.6]; legend = :none;
    sim_plt_1 = plot_fitted_sims_combined(investigation_alt_ps_low_conc.runs[1]; ylimit, yticks, xticks, xguide, lws, ms_true, lw_true, la_true, las, legend)
    sim_plt_2 = plot_fitted_sims_combined(investigation_alt_ps_low_conc.runs[2]; ylimit, yticks, xticks, xguide, lws, ms_true, lw_true, la_true, las, legend)
    sim_plt_3 = plot_fitted_sims_combined(investigation_alt_ps_low_conc.runs[3]; ylimit, yticks, xticks, xguide, lws, ms_true, lw_true, la_true, las, legend)

    # Generates the likelihood profile plots (for each data set).
    no_marker = true; ymin = -4.2; xguide = ""; yguide = ""; yticks = [-4, -2, 0]; xlimit = (-1.0, 1.0); xticks = [-2.0, -1.0, 0.0, 1.0];
    pl_plt_1 = plot_combined_lprofiles(investigation_alt_ps_low_conc.runs[1], :d; lws = [10, 10, 10], no_marker, ymin, xguide, yguide, yticks, las = [0.8, 0.8, 0.8], xlimit, xticks)
    pl_plt_2 = plot_combined_lprofiles(investigation_alt_ps_low_conc.runs[2], :d; lws = [10, 8, 5], no_marker, ymin, xguide, yguide, yticks, las = [0.8, 0.7, 0.6], xlimit, xticks)
    pl_plt_3 = plot_combined_lprofiles(investigation_alt_ps_low_conc.runs[3], :d; lws = [10, 8, 5], no_marker, ymin, xguide, yguide, yticks, las = [0.8, 0.7, 0.6], xlimit, xticks)

    # Generates the ensemble function plots (for each data set).
    xlimit = (0.0, 1.35); ylimit = (0.0, 0.9); xticks = [0.0, 0.45, 0.9, 1.35]; yticks = [0.0, 0.3, 0.6, 0.9]; xguide = ""; lw_true = 13; la_true = 0.6;
    plt_funcs_1 = plot_combined_fitted_funcs(investigation_alt_ps_low_conc.runs[1]; plt_range = (0.0, 1.35), ylimit, xticks, xlimit, yticks, xguide, lw_true, la_true)
    plt_funcs_2 = plot_combined_fitted_funcs(investigation_alt_ps_low_conc.runs[2]; plt_range = (0.0, 1.35), ylimit, xticks, xlimit, yticks, xguide, lw_true, la_true, max_num = 30)
    plt_funcs_3 = plot_combined_fitted_funcs(investigation_alt_ps_low_conc.runs[3]; plt_range = (0.0, 1.35), ylimit, xticks, xlimit, yticks, xguide, lw_true, la_true, max_num = 30)

    # Generates the combined plot.
    bottom_margin = [28mm 28mm 28mm 20mm 20mm 20mm 0mm 0mm 0mm]
    extended_sal_lowconc_ps_plt = plot(sim_plt_1, sim_plt_2, sim_plt_3, pl_plt_1, pl_plt_2, pl_plt_3, plt_funcs_1, plt_funcs_2, plt_funcs_3;
        layout = (3,3), size = (1350, 1200), legendfontsize = 3, left_margin = -5mm, bottom_margin)
    make_fig(extended_sal_lowconc_ps_plt, "Extended_SA_Loop/FigureComponents/extended_sal_pract_ident_lowconc_ps_base")
end

# Main figure - repeat for a different model (GWO).
let
    # Generates the simulation plots (for each data set).
    ylimit = (0.0, 1.2); yticks = [0.0, 0.4, 0.8, 1.2]; xticks = [0.0, 10.0, 20.0, 30.0]; xguide = ""; lws = [10, 8, 6]; ms_true = 7; lw_true = 12; la_true = 1.0; las = [0.8, 0.7, 0.6]; legend = :none;
    sim_plt_1 = plot_fitted_sims_combined(investigation_gwo.runs[1]; ylimit, yticks, xticks, xguide, lws, ms_true, lw_true, la_true, las, legend)
    sim_plt_2 = plot_fitted_sims_combined(investigation_gwo.runs[2]; ylimit, yticks, xticks, xguide, lws, ms_true, lw_true, la_true, las, legend)
    sim_plt_3 = plot_fitted_sims_combined(investigation_gwo.runs[3]; ylimit, yticks, xticks, xguide, lws, ms_true, lw_true, la_true, las, legend)

    # Generates the likelihood profile plots (for each data set).
    no_marker = true; ymin = -4.2; xguide = ""; yguide = ""; yticks = [-4, -2, 0]; xlimit = (-2.0, 0.0); xticks = [-2.0, -1.0, 0.0, 0.0];
    pl_plt_1 = plot_combined_lprofiles(investigation_gwo.runs[1], :dₚ; lws = [10, 10, 10], no_marker, ymin, xguide, yguide, yticks, las = [0.8, 0.8, 0.8], xlimit, xticks)
    pl_plt_2 = plot_combined_lprofiles(investigation_gwo.runs[2], :dₚ; lws = [10, 8, 5], no_marker, ymin, xguide, yguide, yticks, las = [0.8, 0.7, 0.6], xlimit, xticks)
    pl_plt_3 = plot_combined_lprofiles(investigation_gwo.runs[3], :dₚ; lws = [10, 8, 5], no_marker, ymin, xguide, yguide, yticks, las = [0.8, 0.7, 0.6], xlimit, xticks)

    # Generates the ensemble function plots (for each data set).
    xlimit = (0.0, 0.8); ylimit = (0.0, 1.2); xticks = [0.0, 0.2, 0.4, 0.6, 0.8]; yticks = [0.0, 0.4, 0.8, 1.2]; xguide = ""; lw_true = 13; la_true = 0.6;
    plt_funcs_1 = plot_combined_fitted_funcs(investigation_gwo.runs[1]; plt_range = (0.0, 0.8), ylimit, xticks, xlimit, yticks, xguide, lw_true, la_true)
    plt_funcs_2 = plot_combined_fitted_funcs(investigation_gwo.runs[2]; plt_range = (0.0, 0.8), ylimit, xticks, xlimit, yticks, xguide, lw_true, la_true, max_num = 30)
    plt_funcs_3 = plot_combined_fitted_funcs(investigation_gwo.runs[3]; plt_range = (0.0, 0.8), ylimit, xticks, xlimit, yticks, xguide, lw_true, la_true, max_num = 30)

    # Generates the combined plot.
    bottom_margin = [28mm 28mm 28mm 20mm 20mm 20mm 0mm 0mm 0mm]
    gwo_pract_ident_plt = plot(sim_plt_1, sim_plt_2, sim_plt_3, pl_plt_1, pl_plt_2, pl_plt_3, plt_funcs_1, plt_funcs_2, plt_funcs_3;
        layout = (3,3), size = (1350, 1200), legendfontsize = 3, left_margin = -5mm, bottom_margin)
    make_fig(gwo_pract_ident_plt, "Extended_SA_Loop/FigureComponents/gwo_pract_ident_base")
end

# Extended self-activation loop, full identifiability when measuring X only but with lots of data.
let
    ylimit = (0.0, 6.0); yticks = [0.0, 2.0, 4.0, 6.0]; xticks = [0.0, 15.0, 30.0, 45.0]; xguide = ""; lws = [10, 8, 6]; ms_true = 7; lw_true = 12; la_true = 1.0; las = [0.8, 0.7, 0.6]; legend = :none;
    sim_plt = plot_fitted_sims_combined(investigation_high_data; ylimit, yticks,
        xguide, lws, ms_true = 6, lw_true, la_true, las, xticks)
    xlimit = (-1.0, 1.0); xticks = [-2.0, -1.0, 0.0, 1.0];
    pl_plt = plot_combined_lprofiles(investigation_high_data, :d; lw = 10, no_marker = true,
        ymin = -4.2, xguide = "", yguide = "", yticks = [-4, -2, 0], la = 0.7, xlimit, xticks)
    xlimit = (0.1, 4.6); ylimit = (0.0, 1.2); xticks = [0.1, 1.6, 3.1, 4.6]; yticks = [0.0, 0.4, 0.8, 1.2]; xguide = ""; lw_true = 13; la_true = 0.6;
    plt_funcs = plot_combined_fitted_funcs(investigation_high_data; ylimit, xticks, xlimit, yticks, xguide, lw_true, la_true, max_num = 30)
    extended_sal_maxdata_plt = plot(sim_plt, pl_plt, plt_funcs; left_margin = 10Plots.mm, layout = (1,3),
        size = (1350, 350), legendfontsize = 3)
    make_fig(extended_sal_maxdata_plt, "Extended_SA_Loop/FigureComponents/extended_sal_pract_ident_maxdata_base")
end

### Creates Prediction Ensemble Plots ###

# Declares function for generating the individual plots.
function plot_ensemble_predictions(opt_sols, md; n = length(opt_sols), cond = 1, lw = 6,
        color = [1 2], true_color = [:blue3 :red2], lw_true = 4,linestyle = :dash, la = 0.4, kwargs...)
    plot()
    plot!(investigation_main.runs[1].md.sols_true[cond]; color = true_color, lw = lw_true, label = "")
    for idx = 1:n
        fitted_osol = NeuralCRNsIdentifiability.get_fitted_osol(opt_sols[idx]; cond)
        plot!(fitted_osol; lw, color, linestyle, la, label = "", kwargs...)
    end
    NeuralCRNsIdentifiability.plot_data!(md; lw1 = 0, lw2 = 0, ms = 8, kwargs...)
end

# Creates the plot itself.
begin
    ylimit = (0.0, 6.0); yticks = [0.0, 2.0, 4.0, 6.0]; xticks = [0.0, 15.0, 30.0, 45.0]; xguide = ""; legend = :none;
    data1_known_ens_pred_plt = plot_ensemble_predictions(investigation_main.runs[1].opt_sols[1], investigation_main.runs[1].md; ylimit, yticks, xticks, xguide, legend)
    data1_param_ens_pred_plt = plot_ensemble_predictions(investigation_main.runs[1].opt_sols[2], investigation_main.runs[1].md; ylimit, yticks, xticks, xguide, legend)
    data1_nn_ens_pred_plt = plot_ensemble_predictions(investigation_main.runs[1].opt_sols[3], investigation_main.runs[1].md; ylimit, yticks, xticks, xguide, legend)
    data2_known_ens_pred_plt = plot_ensemble_predictions(investigation_main.runs[2].opt_sols[1], investigation_main.runs[2].md; ylimit, yticks, xticks, xguide, legend)
    data2_param_ens_pred_plt = plot_ensemble_predictions(investigation_main.runs[2].opt_sols[2], investigation_main.runs[2].md; ylimit, yticks, xticks, xguide, legend)
    data2_nn_ens_pred_plt = plot_ensemble_predictions(investigation_main.runs[2].opt_sols[3], investigation_main.runs[2].md; ylimit, yticks, xticks, xguide, legend)
    data3_known_ens_pred_plt = plot_ensemble_predictions(investigation_main.runs[3].opt_sols[1], investigation_main.runs[3].md; ylimit, yticks, xticks, xguide, legend)
    data3_param_ens_pred_plt = plot_ensemble_predictions(investigation_main.runs[3].opt_sols[2], investigation_main.runs[3].md; ylimit, yticks, xticks, xguide, legend)
    data3_nn_ens_pred_plt = plot_ensemble_predictions(investigation_main.runs[3].opt_sols[3], investigation_main.runs[3].md; ylimit, yticks, xticks, xguide, legend)

    # Generates the combined plot.
    extended_sal_main_plt = plot(data1_known_ens_pred_plt, data2_known_ens_pred_plt, data3_known_ens_pred_plt, data1_param_ens_pred_plt, data2_param_ens_pred_plt, data3_param_ens_pred_plt, data1_nn_ens_pred_plt, data2_nn_ens_pred_plt, data3_nn_ens_pred_plt;
        layout = (3,3), size = (1350, 1000), legendfontsize = 3, left_margin = -5mm, bottom_margin = -3mm, xguide = "")
    make_fig(extended_sal_main_plt, "Extended_SA_Loop/FigureComponents/extended_sal_pred_uncert_base")
end
