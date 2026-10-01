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
default(framestyle = :box)

# Fetches the completed scan.
begin
    # Fetches the modified SIR misspecification run.
    run_tag_modified_sir = "modified_sir_misspecification_basic_run"
    investigation_modified_sir = load_investigation(run_tag_modified_sir)

    # Fetches the repeat of the modified SIR misspecification run.
    run_tag_modified_sir_repeat = "modified_sir_misspecification_repeat_run"
    investigation_modified_sir_repeat = load_investigation(run_tag_modified_sir_repeat)

    # Fetches the incoherent feedforward loop misspecification run.
    run_tag_ic_feedforward = "ic_feedforward_misspecification_basic_run"
    investigation_ic_feedforward = load_investigation(run_tag_ic_feedforward)
end


# Fetches the completed runs.
@time begin
    # Fetches the modified SIR misspecification run.
    run_sir_1_1 = "modified_sir_misspecification_basic_run_1_1"
    investigation_sir_1_1 = load_investigation(run_sir_1_1)
    run_sir_1_2 = "modified_sir_misspecification_basic_run_1_2"
    investigation_sir_1_2 = load_investigation(run_sir_1_2)
    run_sir_2_1 = "modified_sir_misspecification_basic_run_2_1"
    investigation_sir_2_1 = load_investigation(run_sir_2_1)
    run_sir_2_2 = "modified_sir_misspecification_basic_run_2_2"
    investigation_sir_2_2 = load_investigation(run_sir_2_2)

    # Fetches the incoherent feedforward loop misspecification run.
    run_ic_1_1 = "ic_feedforward_misspecification_basic_run_1_1"
    investigation_ic_1_1 = load_investigation(run_ic_1_1)
    run_ic_1_2 = "ic_feedforward_misspecification_basic_run_1_2"
    investigation_ic_1_2 = load_investigation(run_ic_1_2)
    run_ic_2_1 = "ic_feedforward_misspecification_basic_run_2_1"
    investigation_ic_2_1 = load_investigation(run_ic_2_1)
    run_ic_2_2 = "ic_feedforward_misspecification_basic_run_2_2"
    investigation_ic_2_2 = load_investigation(run_ic_2_2)

    # Fetches other runs.
    run_ic_basic = "ic_feedforward_misspecification_basic_run"
    investigation_ic_basic = load_investigation(run_ic_basic)
    run_ic_short = "ic_feedforward_misspecification_basic_run_short"
    investigation_ic_short = load_investigation(run_ic_short)
    run_sir_basic = "modified_sir_misspecification_basic_run"
    investigation_sir_basic = load_investigation(run_sir_basic)
    run_sir_short = "modified_sir_misspecification_basic_run_short"
    investigation_sir_short = load_investigation(run_sir_short)
    run_sir_extra = "modified_sir_misspecification_extra_run"
    investigation_sir_extra = load_investigation(run_sir_extra)
    run_sir_extra_short = "modified_sir_misspecification_extra_run_short"
    investigation_sir_extra_short = load_investigation(run_sir_extra_short)
end
# Fetches the completed runs.
@time begin
    # Fetches the modified SIR misspecification run.
    run_sir_r1 = "modified_sir_misspecification_basic_run_2_1_r1"
    investigation_sir_r1 = load_investigation(run_sir_r1)
    run_sir_r2 = "modified_sir_misspecification_basic_run_2_1_r2"
    investigation_sir_r2 = load_investigation(run_sir_r2)
    run_sir_r3 = "modified_sir_misspecification_basic_run_2_1_r3"
    investigation_sir_r3 = load_investigation(run_sir_r3)
end

### Figure generation ###

# Functions for generating the figures.
function make_sir_misspec_plot(run)
    legend = :none
    bottom_margin = [-2mm -2mm -2mm -2mm 0mm -2mm]
    left_margin = [-3mm -3mm -3mm -3mm 9mm 9mm]
    return plot(
        plot_all_fitted_simulation(run, 1; xlimit = (0.0, 1000.0), xticks = [0.0, 250.0, 500.0, 750.0, 1000.0], ylimit = (0.0, 100.0), yticks = [0.0, 25, 50, 75, 100.0]),
        plot_all_fitted_simulation(run, 2; xlimit = (0.0, 1000.0), xticks = [0.0, 250.0, 500.0, 750.0, 1000.0], ylimit = (0.0, 100.0), yticks = [0.0, 25, 50, 75, 100.0]),
        plot_all_fitted_simulation(run, 3; xlimit = (0.0, 1000.0), xticks = [0.0, 250.0, 500.0, 750.0, 1000.0], ylimit = (0.0, 100.0), yticks = [0.0, 25, 50, 75, 100.0]),
        plot_all_fitted_simulation(run, 4; xlimit = (0.0, 1000.0), xticks = [0.0, 250.0, 500.0, 750.0, 1000.0], ylimit = (0.0, 100.0), yticks = [0.0, 25, 50, 75, 100.0]),
        plot(; legend=false, grid=false, ticks=nothing, framestyle=:none, background_color=:transparent),
        plot_combined_lprofiles(run, :γ; no_marker = true, ymin = -4.2, xguide = "", yguide = "", yticks = [-4, -2, 0], xlimit = (-2.5, -2.2), xticks = [-2.5, -2.4, -2.3, -2.2]);
        bottom_margin, left_margin, legend, xguide = "", yguide = "",
        size = (1200, 600), layout = @layout [grid(2,2){0.65w} [a; b{0.6h}]]
    )
end
function make_ic_misspec_plot(run)
    legend = :none
    bottom_margin = [-2mm -2mm -2mm -2mm 0mm -2mm]
    left_margin = [-3mm -3mm -3mm -3mm 9mm 9mm]
    return plot(
        plot_all_fitted_simulation(run, 1; xlimit = (0.0, 10.0), xticks = [0.0, 2.5, 5.0, 7.5, 10.0], ylimit = (0.0, 3.2), yticks = [0.0, 0.8, 1.6, 2.4, 3.2]),
        plot_all_fitted_simulation(run, 2; xlimit = (0.0, 10.0), xticks = [0.0, 2.5, 5.0, 7.5, 10.0], ylimit = (0.0, 3.2), yticks = [0.0, 0.8, 1.6, 2.4, 3.2]),
        plot_all_fitted_simulation(run, 3; xlimit = (0.0, 10.0), xticks = [0.0, 2.5, 5.0, 7.5, 10.0], ylimit = (0.0, 3.2), yticks = [0.0, 0.8, 1.6, 2.4, 3.2]),
        plot_all_fitted_simulation(run, 4; xlimit = (0.0, 10.0), xticks = [0.0, 2.5, 5.0, 7.5, 10.0], ylimit = (0.0, 3.2), yticks = [0.0, 0.8, 1.6, 2.4, 3.2]),
        plot(; legend=false, grid=false, ticks=nothing, framestyle = :none, background_color=:transparent),
        plot_combined_lprofiles(run, :d; no_marker = true, colors = [:darkorange :purple :green :cyan], ymin = -4.2, xguide = "", yguide = "", yticks = [-4, -2, 0], xlimit = (-1.5, 0.0), xticks = [-1.5, -1.0, -0.5, 0.0]);
        bottom_margin, left_margin, legend, xguide = "", yguide = "",
        size = (1200, 600), layout = @layout [grid(2,2){0.65w} [a; b{0.6h}]]
    )
end

# Generates the main figure.
modified_sir_plt = make_sir_misspec_plot(investigation_modified_sir)
make_fig(modified_sir_plt, "Modified_SIR/FigureComponents/modified_sir_misspecification_base")

# Generates the repeat figure.
modified_sir_plt_repeat = make_sir_misspec_plot(investigation_modified_sir_repeat)
make_fig(modified_sir_plt_repeat, "Modified_SIR/FigureComponents/modified_sir_misspecification_repeat_base")

# Generates the ic model figure.
ic_feedforward_plt = make_ic_misspec_plot(investigation_ic_feedforward)
make_fig(ic_feedforward_plt, "Modified_SIR/FigureComponents/ic_feedforward_misspecification_base")


# Generates the main figure.
let
    legend = :none
    bottom_margin = [-2mm -2mm -2mm -2mm 0mm -2mm]
    left_margin = [-3mm -3mm -3mm -3mm 9mm 9mm]
    modified_sir_plt = plot(
        plot_all_fitted_simulation(investigation_modified_sir, 1; xlimit = (0.0, 1000.0), xticks = [0.0, 250.0, 500.0, 750.0, 1000.0], ylimit = (0.0, 100.0), yticks = [0.0, 25, 50, 75, 100.0]),
        plot_all_fitted_simulation(investigation_modified_sir, 2; xlimit = (0.0, 1000.0), xticks = [0.0, 250.0, 500.0, 750.0, 1000.0], ylimit = (0.0, 100.0), yticks = [0.0, 25, 50, 75, 100.0]),
        plot_all_fitted_simulation(investigation_modified_sir, 3; xlimit = (0.0, 1000.0), xticks = [0.0, 250.0, 500.0, 750.0, 1000.0], ylimit = (0.0, 100.0), yticks = [0.0, 25, 50, 75, 100.0]),
        plot_all_fitted_simulation(investigation_modified_sir, 4; xlimit = (0.0, 1000.0), xticks = [0.0, 250.0, 500.0, 750.0, 1000.0], ylimit = (0.0, 100.0), yticks = [0.0, 25, 50, 75, 100.0]),
        plot(; legend=false, grid=false, ticks=nothing, framestyle=:none, background_color=:transparent),
        plot_combined_lprofiles(investigation_modified_sir, :γ; no_marker = true, ymin = -4.2, xguide = "", yguide = "", yticks = [-4, -2, 0], xlimit = (-2.5, -2.2), xticks = [-2.5, -2.4, -2.3, -2.2]);
        bottom_margin, left_margin, legend, xguide = "", yguide = "",
        size = (1200, 600), layout = @layout [grid(2,2){0.65w} [a; b{0.6h}]]
    )
    make_fig(modified_sir_plt, "Modified_SIR/FigureComponents/modified_sir_misspecification_base")
end

# Generates the IC feedforward loop supplementary figure.
let
    legend = :none
    bottom_margin = [-2mm -2mm -2mm -2mm 0mm -2mm]
    left_margin = [-3mm -3mm -3mm -3mm 9mm 9mm]
    ic_feedforward_plt = plot(
        plot_all_fitted_simulation(investigation_ic_feedforward, 1; xlimit = (0.0, 10.0), xticks = [0.0, 2.5, 5.0, 7.5, 10.0], ylimit = (0.0, 3.2), yticks = [0.0, 0.8, 1.6, 2.4, 3.2]),
        plot_all_fitted_simulation(investigation_ic_feedforward, 2; xlimit = (0.0, 10.0), xticks = [0.0, 2.5, 5.0, 7.5, 10.0], ylimit = (0.0, 3.2), yticks = [0.0, 0.8, 1.6, 2.4, 3.2]),
        plot_all_fitted_simulation(investigation_ic_feedforward, 3; xlimit = (0.0, 10.0), xticks = [0.0, 2.5, 5.0, 7.5, 10.0], ylimit = (0.0, 3.2), yticks = [0.0, 0.8, 1.6, 2.4, 3.2]),
        plot_all_fitted_simulation(investigation_ic_feedforward, 4; xlimit = (0.0, 10.0), xticks = [0.0, 2.5, 5.0, 7.5, 10.0], ylimit = (0.0, 3.2), yticks = [0.0, 0.8, 1.6, 2.4, 3.2]),
        plot(; legend=false, grid=false, ticks=nothing, framestyle = :none, background_color=:transparent),
        plot_combined_lprofiles(investigation_ic_feedforward, :d; no_marker = true, colors = [:darkorange :purple :green :cyan], ymin = -4.2, xguide = "", yguide = "", yticks = [-4, -2, 0], xlimit = (-1.5, 0.0), xticks = [-1.5, -1.0, -0.5, 0.0]);
        bottom_margin, left_margin, legend, xguide = "", yguide = "",
        size = (1200, 600), layout = @layout [grid(2,2){0.65w} [a; b{0.6h}]]
    )
    make_fig(ic_feedforward_plt, "Modified_SIR/FigureComponents/ic_feedforward_misspecification_base")
end
