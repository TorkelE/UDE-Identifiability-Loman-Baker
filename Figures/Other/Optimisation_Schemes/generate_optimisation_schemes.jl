### Preparations ###

# Sets the correct directory. Fetches packages.
begin
    cd("$(@__DIR__)/../../..")
    using Pkg; Pkg.activate(".")
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

# Fetch various sets of optimisers.
begin
    using OptimizationNLopt
    using OptimizationOptimJL
    using OptimizationOptimisers: Adam
end

# Set plotting defaults. Creates the (Stable) RNG.
begin
    default(grid = false, framestyle = :box)
    rng = StableRNG(1)
end

### Optimisation Preparation ###


# Sets the true parameter values (and their ranges).
# Beware to make sure log/not-log scale is correct for each parameter.
# Loads the model (we will only use the UDE one for optimisation, the known one for data generation).
begin
    ps_true = [
        :v => log10(1.1),
        :K => log10(2.0),
        :n => 3.0,
        :d => log10(0.5)
    ]
    ps_ranges, (rs_known, _, rs_nn) = get_extended_sal(ps_true, rng)
end


# Generates the synthetic data and plots it.
begin
    u0s = [
        [:X => 2.0, :Y => 0.1]
    ]
    md = MeasuredData(
        rs_known, #The base model.
        u0s, # The initial conditions.
        45.0, # The simulation time.
        ps_true, # The true parameter set.
        20, # The number of samples from each simulation.
        [:X, :Y], # The measured species.
        0.1; # The noise level.
        rng, # The random number generator.
        ps_ranges # The valid ranges assigned to each parameter.
    )
    plot_all_data(md)
end

# Creates an `OptimizationProblem` for the UDE model.
begin
    @everywhere n = 40
    @everywhere opt_prob = $(make_eoptprob(rs_nn, md, n; rng, success_required = false, max_loss_thres = 1e5, min_loss_thres = -Inf))
end

# Fits the model using various optimisers.
begin
    optalg = Adam(5e-3)
    @time opt_sol_adam = fit_model(opt_prob, n; optalg, ealg = EnsembleDistributed(), maxiters = 10000)
    summarise_optimisation(opt_sol_adam, md)
end
begin
    optalg = OptimizationOptimJL.LBFGS()
    @time opt_sol_lbfgs = fit_model(opt_prob, n; optalg, ealg = EnsembleDistributed(), maxiters = 1)
    summarise_optimisation(opt_sol_lbfgs, md)
end
begin
    optalg = NLopt.LN_NELDERMEAD()
    @time opt_sol_neldermead = fit_model(opt_prob, n; optalg, ealg = EnsembleDistributed(), maxiters = 10000)
    summarise_optimisation(opt_sol_neldermead, md)
end
begin
    optalg = Optim.ConjugateGradient()
    @time opt_sol_conjgrad = fit_model(opt_prob, n; optalg, ealg = EnsembleDistributed(), maxiters = 10000)
    summarise_optimisation(opt_sol_conjgrad, md)
end

### Generates Plots ###


# Helper function for plotting loss traces.
function specialised_plot_loss_traces(opt_sol, md, best_loss; include_tuning = false, la = 0.5, lw = 4, is_lbfgs = false, kwargs...)
    plot()
    sols = filter(sol -> sol.opt_sol1.retcode != ReturnCode.Failure, opt_sol.u)
    is_lbfgs && (sols = opt_sol.u)
    for sol in sols
        if !is_lbfgs
            its_main = length(sol.min_trace) - sol.opt_sol2.stats.iterations - 1
            plot!(1:its_main, 1 .- best_loss .+ sol.min_trace[1:its_main]; label = "", lw, la, color = 1)
            if include_tuning
                its_tuning = length(sol.min_trace) - its_main
                plot!(its_main .+ (1:its_tuning), 1 .- best_loss .+ sol.min_trace[its_main .+ (1:its_tuning)]; label = "", lw, la, color = 2)
            end
        else
            plot!(1:length(sol.min_trace[2:end]), 1 .- best_loss .+ sol.min_trace[2:end]; label = "", lw, la, color = 1)
        end
    end
    plot!(yaxis = :log10, ylimit = (1.0, 1e4), xlimit = (1, 11000))
end


begin # First run this.
    ylimit = (0.0, 6.0); yticks = [0.0, 2.0, 4.0, 6.0]; xticks = [0.0, 15.0, 30.0, 45.0]; xguide = ""; ms = 5; legend = :none;
    sim_plt_1 = plot_fitted_simulation(opt_sol_adam, md; ylimit, yticks, xticks, xguide, lws, ms, legend)
    sim_plt_2 = plot_fitted_simulation(opt_sol_lbfgs, md; ylimit, yticks, xticks, xguide, lws, ms, legend)
end
begin # While using Revise, change `get_fitted_oprob` to use `best_u1`.
    sim_plt_3 = plot_fitted_simulation(opt_sol_neldermead, md; ylimit, yticks, xticks, xguide, lws, ms, legend)
    sim_plt_4 = plot_fitted_simulation(opt_sol_conjgrad, md; ylimit, yticks, xticks, xguide, lws, ms, legend)
end
begin # While using Revise, change `get_fitted_oprob` back to use `best_u`.
    yticks = [10e0, 10e1, 10e2, 10e3, 10e4]; ylimit = (10e0, 10e4); xticks = [0.0, 2500, 5000, 7500, 10000];
    opt_plt_1 = specialised_plot_loss_traces(opt_sol_adam, md, best_loss; include_tuning = true, xticks, ylimit, yticks)
    opt_plt_2 = specialised_plot_loss_traces(opt_sol_lbfgs, md, best_loss; is_lbfgs = true, xticks, ylimit, yticks)
    opt_plt_3 = specialised_plot_loss_traces(opt_sol_neldermead, md, best_loss; xticks, ylimit, yticks)
    opt_plt_4 = specialised_plot_loss_traces(opt_sol_conjgrad, md, best_loss; xticks, ylimit, yticks)
    opt_choice_plot = plot(
        sim_plt_1, sim_plt_2, sim_plt_3, sim_plt_4,
        opt_plt_1, opt_plt_2, opt_plt_3, opt_plt_4;
        layout = (2, 4),
        size = (1600, 650),
        legend = false, bottom_margin = 8mm, left_margin = -4mm
    )
    make_fig(opt_choice_plot, "Other/Optimisation_Schemes/FigureComponents/optimiser_choice_base")
end
