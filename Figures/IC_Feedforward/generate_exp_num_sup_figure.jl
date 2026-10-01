### Preparations ###

# Sets the correct directory. Fetches packages.
begin
    cd("$(@__DIR__)/../..")
    using Pkg; Pkg.activate(".")
    using NeuralCRNsIdentifiability
    using Plots
    import Plots.mm
end

# Sets up multi-process parallelisation.
begin
    using Distributed
    addprocs(42)
end

# Fetch packages.
begin
    @everywhere using NeuralCRNsIdentifiability
    using Plots
    using Optimization
    using StableRNGs
    import Plots.mm
end

# Set plotting defaults. Creates the (Stable) RNG.
begin
    default(grid = false, framestyle = :box)
    rng = StableRNG(12345678)
end

# Sets the true parameter values (and their ranges).
# Beware to make sure log/not-log scale is correct for each parameter.
# Loads the models and the parameter ranges.
begin
    ps_true = [
        :d => log10(0.5),
        :v1 => log10(5.0),
        :K1 => log10(0.8),
        :n => 5,
        :v2 => log10(1.0),
        :K2 => log10(0.6)
        ]
    ps_ranges, (rs_known, rs_params_params, rs_nn_params, rs_params_nn, rs_nn_nn) = get_ic_feedforward(ps_true, rng)
end

### Data Generation ###

# Generates the synthetic datasets and plots them
begin
    u0s = [
        [:X => 0.0, :Y => 0.0, :Z => 0.0],
        [:X => 1.0, :Y => 0.0, :Z => 0.0],
        [:X => 3.0, :Y => 3.0, :Z => 0.0],
        [:X => 0.0, :Y => 3.0, :Z => 0.0],
    ]
    md1 = MeasuredData(
        rs_known, #The base model.
        u0s, # The initial conditions.
        10.0, # The simulation time.
        ps_true, # The true parameter set.
        30, # The number of samples from each simulation.
        [:Z], # The measured species.
        0.15; # The noise level.
        rng, # The random number generator.
        ps_ranges # The valid ranges assigned to each parameter.
    )
    plot_all_data(md1)
end
begin
    u0s = [
        [:X => 0.0, :Y => 0.0, :Z => 0.0]
    ]
    md2 = MeasuredData(
        rs_known, #The base model.
        u0s, # The initial conditions.
        10.0, # The simulation time.
        ps_true, # The true parameter set.
        30, # The number of samples from each simulation.
        [:Z], # The measured species.
        0.15; # The noise level.
        rng, # The random number generator.
        ps_ranges # The valid ranges assigned to each parameter.
    )
    plot_all_data(md2)
end
begin
    u0s = [
        [:X => 0.0, :Y => 0.0, :Z => 0.0]
    ]
    md3 = MeasuredData(
        rs_known, #The base model.
        u0s, # The initial conditions.
        10.0, # The simulation time.
        ps_true, # The true parameter set.
        30, # The number of samples from each simulation.
        [:Y, :Z], # The measured species.
        0.15; # The noise level.
        rng, # The random number generator.
        ps_ranges # The valid ranges assigned to each parameter.
    )
    plot_all_data(md3)
end

# Creates an solves an ensemble optimization problems for n multistarts.
begin
    rs = rs_nn_nn
    @everywhere n = 40
    @everywhere opt_prob1 = $(make_eoptprob(rs, md1, n; rng))
    @time opt_sol1 = fit_model(opt_prob1, n; ealg = EnsembleDistributed(), maxiters = 2000)
    summarise_optimisation(opt_sol1, md1)
end
begin
    rs = rs_nn_nn
    @everywhere n = 40
    @everywhere opt_prob2 = $(make_eoptprob(rs, md2, n; rng))
    @time opt_sol2 = fit_model(opt_prob2, n; ealg = EnsembleDistributed(), maxiters = 2000)
    summarise_optimisation(opt_sol2, md2)
end
begin
    rs = rs_nn_nn
    @everywhere n = 40
    @everywhere opt_prob3 = $(make_eoptprob(rs, md3, n; rng))
    @time opt_sol3 = fit_model(opt_prob3, n; ealg = EnsembleDistributed(), maxiters = 2000)
    summarise_optimisation(opt_sol3, md3)
end

### Generate Plot ###
begin
    xlimit = (0.0, 10.0); xticks = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0]; ylimit = (0.0, 5.0); yticks = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
    p1 = plot_all_fitted_simulation(opt_sol1, md1; xlimit, xticks, ylimit, yticks)
    p2 = plot_all_fitted_simulation(opt_sol2, md2; xlimit, xticks, ylimit, yticks) # Color order is weird (due to routine in scripts where it was easier to simply recolor this one case manually to use the same colour scheme as in the other two figures).
    p3 = plot_all_fitted_simulation(opt_sol3, md3; xlimit, xticks, ylimit, yticks)
    ic_exp_num_analysis = plot(
        p2, p3, p1;
        size = (1320, 700),
        bottom_margin = -4mm,
        layout = @layout [[a; b] c{0.66w}]
    )
    make_fig(ic_exp_num_analysis, "IC_Feedforward/FigureComponents/ic_feedforward_exp_num_analysis_base")
end
