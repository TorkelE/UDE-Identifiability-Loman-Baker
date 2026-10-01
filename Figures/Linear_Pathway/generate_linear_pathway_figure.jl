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
    # Linear pathway scan.
    run_tag_ln_pathway = "linear_pathway_nn_constraints_scan"
    check_scan(run_tag_ln_pathway; print_info = true)
    investigation_ln_pathway = load_investigation(run_tag_ln_pathway)

    # Goodwind oscillator scan.
    run_tag_gw_oscillator = "gwo_nn_constraints_scan"
    check_scan(run_tag_gw_oscillator; print_info = true)
    investigation_gw_oscillator = load_investigation(run_tag_gw_oscillator)
end

### Figure generation ###

# Plotting function.
function make_nn_constraints_plot(run; ms = 5, legend = :none, plt_range = (0.0, 1.2), model = :ln, la_true = 0.5, kwargs...)
    # Deterimnes plotting settings.
    if model == :ln
        p = :d
        xlimit1 = (0.0, 5.0); xlimit2 = (-0.5, 1.0); xlimit3 = (0.0, 1.2);
        ylimit1 = (0.0, 1.5);  ylimit3 = (0.0, 1.6);
        xticks1 = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]; xticks2 = [-0.5, 0.0, 0.5, 1.0]; xticks3 = [0.0, 0.6, 1.2];
        yticks1 = [0.0, 0.5, 1.0, 1.5]; yticks2 = [-4, -2, 0]; yticks3 = [0.0, 0.4, 0.8, 1.2, 1.6];
    elseif model == :gwo
        p = :dₚ
        xlimit1 = (0.0, 30.0); xlimit2 = (-2.0, -0.5); xlimit3 = (0.0, 0.8);
        ylimit1 = (0.0, 1.5);  ylimit3 = (0.0, 1.6);
        xticks1 = [0.0, 7.5, 15.0, 22.5, 30.0]; xticks2 = [-2, -1.5, -1.0, -0.5]; xticks3 = [0.0, 0.4, 0.8];
        yticks1 = [0.0, 0.5, 1.0, 1.5]; yticks2 = [-4, -2, 0]; yticks3 = [0.0, 0.4, 0.8, 1.2, 1.6]
    else
        error("Unknown model.")
    end

    # Make individual plots.
    plt1 = plot_all_fitted_simulation(run, 2; ms, ylimit = ylimit1, yticks = yticks1, xlimit = xlimit1, xticks = xticks1)
    plt2 = plot_all_fitted_simulation(run, 3; ms, ylimit = ylimit1, yticks = yticks1, xlimit = xlimit1, xticks = xticks1)
    plt3 = plot_all_fitted_simulation(run, 4; ms, ylimit = ylimit1, yticks = yticks1, xlimit = xlimit1, xticks = xticks1)
    plt4 = plot_all_fitted_simulation(run, 5; ms, ylimit = ylimit1, yticks = yticks1, xlimit = xlimit1, xticks = xticks1)
    plt5 = plot_combined_lprofiles(run, p; ymin = -4.2, xguide = "", yguide = "", yticks = yticks2, xlimit = xlimit2, xticks = xticks2)
    plt6 = plot([], []; xlabel = "", ylabel = "", title = "", legend = false, framestyle = :none)
    plt7 = plot_fitted_funcs(run, 2; plt_range, ylimit = ylimit3, yticks = yticks3, xlimit = xlimit3, xticks = xticks3)
    plt8 = plot_fitted_funcs(run, 3; plt_range, ylimit = ylimit3, yticks = yticks3, xlimit = xlimit3, xticks = xticks3)
    plt9 = plot_fitted_funcs(run, 4; plt_range, ylimit = ylimit3, yticks = yticks3, xlimit = xlimit3, xticks = xticks3, la_true)
    plt10 = plot_fitted_funcs(run, 5; plt_range, ylimit = ylimit3, yticks = yticks3, xlimit = xlimit3, xticks = xticks3)

    # Composes the final plot.
    bottom_margin = [-2mm -2mm 6mm 6mm 0mm 0mm 0mm 0mm 0mm 0mm 0mm]
    left_margin = [0mm 0mm -4mm -4mm 12mm 12mm 12mm 0mm -4mm -4mm -4mm]
    plot(
        plt1, plt2, plt3, plt4, plt6, plt5, plt6, plt7, plt8, plt9, plt10;
        legend, left_margin, bottom_margin, xguide = "", yguide = "",
        size = (1200, 830), layout = @layout [[grid(2,2){0.6w} [a{0.26h}; b{0.62h}; c]]; grid(1,4){0.3h}]
    )
end

# Generates the figures.
ln_plt1 = make_nn_constraints_plot(investigation_ln_pathway.runs[1]; la_true = 0.75)
ln_plt2 = make_nn_constraints_plot(investigation_ln_pathway.runs[2]; la_true = 0.75)
ln_plt3 = make_nn_constraints_plot(investigation_ln_pathway.runs[3]; la_true = 0.75)
gwo_plt1 = make_nn_constraints_plot(investigation_gw_oscillator.runs[1]; model = :gwo, la_true = 0.75)
gwo_plt2 = make_nn_constraints_plot(investigation_gw_oscillator.runs[2]; model = :gwo, la_true = 0.75)

# Saves the figures.
make_fig(ln_plt1, "Linear_Pathway/FigureComponents/ln_pathway_nn_constraints_repeat1_base")
make_fig(ln_plt2, "Linear_Pathway/FigureComponents/ln_pathway_nn_constraints_repeat2_base")
make_fig(ln_plt3, "Linear_Pathway/FigureComponents/ln_pathway_nn_constraints_repeat3_base")
make_fig(gwo_plt1, "Linear_Pathway/FigureComponents/gwo_nn_constraints_repeat1_base")
make_fig(gwo_plt2, "Linear_Pathway/FigureComponents/gwo_nn_constraints_repeat2_base")

