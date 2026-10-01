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

# Fetches the completed scan.
begin
    run_tag = "simple_sal_basic_run"
    loaded_investigation = load_investigation(run_tag)
end

### Figure generation ###

# Generates and saves the simple self-activation loop figure.
let
    sim_plt = plot_fitted_sims_combined(loaded_investigation; ylimit = (0.0, 1.3), yticks = [0.0, 0.6, 1.2],
        xguide = "", lws = [10,8,6], ms_true = 7, lw_true = 12, la_true = 1.0, las = [0.8, 0.7, 0.6])
    pl_plt = plot_combined_lprofiles(loaded_investigation, :d; lw = 10, no_marker = true,
        ymin = -4.2, xguide = "", yguide = "", yticks = [-4, -2, 0], la = 0.7)
    plt_funcs = plot_combined_fitted_funcs(loaded_investigation; ylimit = (0.0, 1.0),
        xlimit = (0.15, 1.35), xticks = [0.15, 0.55, 0.95, 1.35], yticks = [0.0, 0.25, 0.5, 0.75, 1.0],
        xguide = "", lw_true = 13, la_true = 0.6)
    plt_sa_loop = plot(sim_plt, pl_plt, plt_funcs; left_margin = 10mm, layout = (1,3),
        size = (1350, 350), legendfontsize = 3)
    make_fig(plt_sa_loop, "Simple_SA_Loop/FigureComponents/simple_sal_struc_nonident_base")
end
