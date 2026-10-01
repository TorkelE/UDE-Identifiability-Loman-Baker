### Preparations ###

# Sets the correct directory. Fetches packages.
begin
    cd("$(@__DIR__)/../..")
    using Pkg; Pkg.activate(".")
    using NeuralCRNsIdentifiability
    using Plots
    import Plots.mm
end

# Set plotting defaults.
default(framestyle = :box, grid = true)

# Fetches the completed scan.
begin
    run_tag_ic_feedforward = "ic_feedforward_num_nn_variation_basic_scan"
    check_scan(run_tag_ic_feedforward; print_info = true)
    investigation_ic_feedforward = load_investigation(run_tag_ic_feedforward)
end


### Figure generation ###

# Function for generating the plot for a specific run.
function plot_ic_run(run; legend = :none, xlimit1 = (0.0, 2.25), ylimit1 = (0.0, 6.0), xticks1 = [0.0, 0.75, 1.5, 2.25], yticks1 = [0.0, 2.0, 4.0, 6.0],
        xlimit2 = (0.0, 2.25), ylimit2 = (0.0, 1.2), xticks2 = [0.0, 0.75, 1.6, 2.25], yticks2 = [0.0, 0.4, 0.8, 1.2], ms_true = 7)
    bottom_margin = [45mm 45mm -6mm -6mm 0mm 0mm -6mm -6mm 0mm 0mm]
    left_margin = [0mm 10mm 0mm -7mm 0mm -7mm 8mm -7mm 8mm -7mm]
    plot(
        plot_fitted_sims_combined(run; legend = :none, ms_true, xlimit = (0.0, 8.0), xticks = [0.0, 2.0, 4.0, 6.0, 8.0], ylimit = (0.0, 2.7), yticks = [0.0, 0.9, 1.8, 2.7], plot_models = 2:5, linestyles = [:dot, :dash, :dashdotdot, :dashdot], colors = plotting_colors_ic),
        plot_combined_lprofiles(run, :d; no_marker = true, ymin = -4.2, xguide = "", yguide = "", yticks = [-4, -2, 0], xlimit = (-1.0, 1.0), xticks = [-2.0, -1.0, 0.0, 1.0]),
        plot(
            plot(
                plot_fitted_funcs(run, 2; f_idx = 1, xlimit = xlimit1, ylimit = ylimit1, xticks = xticks1, yticks = yticks1, legend),
                plot_fitted_funcs(run, 3; f_idx = 1, xlimit = xlimit1, ylimit = ylimit1, xticks = xticks1, yticks = yticks1, legend),
                plot_fitted_funcs(run, 4; f_idx = 1, xlimit = xlimit1, ylimit = ylimit1, xticks = xticks1, yticks = yticks1, legend),
                plot_fitted_funcs(run, 5; f_idx = 1, xlimit = xlimit1, ylimit = ylimit1, xticks = xticks1, yticks = yticks1, legend);
                layout = (2,2), size = (660, 400)
            ),
            plot(
                plot_fitted_funcs(run, 2; f_idx = 2, xlimit = xlimit2, ylimit = ylimit2, xticks = xticks2, yticks = yticks2, legend),
                plot_fitted_funcs(run, 3; f_idx = 2, xlimit = xlimit2, ylimit = ylimit2, xticks = xticks2, yticks = yticks2, legend),
                plot_fitted_funcs(run, 4; f_idx = 2, xlimit = xlimit2, ylimit = ylimit2, xticks = xticks2, yticks = yticks2, legend),
                plot_fitted_funcs(run, 5; f_idx = 2, xlimit = xlimit2, ylimit = ylimit2, xticks = xticks2, yticks = yticks2, legend);
                layout = (2,2), size = (660, 400)
            );
            layout = (1,2), size = (1320, 400)
        );
        xguide = "", yguide = "", bottom_margin, left_margin, size = (1320, 1110), layout = @layout [[a{0.6w} b]; c{0.5h}]
    )
end

# Generates the main figure.
ic_feedforward_plt = plot_ic_run(investigation_ic_feedforward.runs[1])
make_fig(ic_feedforward_plt, "IC_Feedforward/FigureComponents/ic_feedforward_2nn_base")

# Generates the repeat of the main figure.
ic_feedforward_repeat_plt = plot_ic_run(investigation_ic_feedforward.runs[2])
make_fig(ic_feedforward_repeat_plt, "IC_Feedforward/FigureComponents/ic_feedforward_2nn_repeat_base")

