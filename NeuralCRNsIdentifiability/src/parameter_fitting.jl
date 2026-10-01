# Functions for fitting the parameters of a model to data.

### Optimisation-Related Structures ###

# An expanded optimization solution, containing additional information.
struct ExpandedOptimizationSolution{S,T}
    """The final, optimal, solution. When monotonicity have been enforced, this has been transformed."""
    u::Vector{Float64}
    """The objective value of the final solution."""
    objective::Float64
    """The optimisation solution of the initial (Adam) optimisation run."""
    opt_sol1::S
    """The optimisation solution of the second (LBFGS) optimisation run."""
    opt_sol2::T
    """The loss function trace of the full optimisation process."""
    trace::Vector{Float64}
    """The minimum loss function trace so far of the full optimisation process."""
    min_trace::Vector{Float64}

    function ExpandedOptimizationSolution(opt_sol1::S, opt_sol2::T, trace::Vector{Float64};
            final_prints = true) where {S,T}
        # Performs final checks on the solution.
        final_prints && (1.1*opt_sol1.objective < opt_sol2.objective) && @warn "The LBFGS tuning did not improve the solution ($(opt_sol1.objective) vs $(opt_sol2.objective))"
        final_prints && println("Adam loss: $(opt_sol1.objective) \t LBFGS loss: $(opt_sol2.objective)")

        # Finds the optimal solution, and saves the corresponding solution and objective value.
        if opt_sol1.objective < opt_sol2.objective
            u = opt_sol1.u
            objective = opt_sol1.objective
        else
            u = opt_sol2.u
            objective = opt_sol2.objective
        end
        if opt_sol1.prob.p[10] != :false
            u = enforce_nn_monotonicity(u, opt_sol1.prob.p[10], opt_sol1.prob.p[11])
        end

        # Finds the trace, returns the updated optimisation solution.
        min_trace = zeros(length(trace)); min_trace[1] = trace[1];
        foreach(i -> min_trace[i] = min(min_trace[i-1], trace[i]), 2:length(min_trace))
        new{S,T}(u, objective, opt_sol1, opt_sol2, trace, min_trace)
    end
end

### Optimisation Problem Generation ###

# Generates an ensemble `OptimizationProblem` for multistart optimization. Handles
# - Selection of initial parameter value guesses (and resamples of bad guesses).
# - Set-up of the callback.
# - The saving of the optimisation traces.
# - Set-up of the "finalisinG" LBFGS run for each optimisation run.
function make_eoptprob(rs::ReactionSystem, md::MeasuredData, n::Int64; rng = Random.default_rng(),
        pinits = make_pinits(rs, md, n, rng), final_prints = true, cb_prog_print = false,
        max_loss_thres = 1e10, success_required = false, rerun_print = true, min_loss_thres = true_loss(md))
    opt_prob_base = make_optprob(rs, md; rng)
    traces = [Float64[] for _ in 1:n]
    prob_func = make_prob_func_eopt(rs, md, n, traces, opt_prob_base; cb_prog_print, max_loss_thres, pinits, min_loss_thres)
    output_func = make_output_func_eopt(traces, opt_prob_base, md; final_prints, max_loss_thres, success_required, rerun_print)
    return EnsembleProblem(opt_prob_base; prob_func, output_func, safetycopy = false)
end

# From a model and some measured data, creates a basic optimisation problem.
function make_optprob(rs::ReactionSystem, md::MeasuredData;
        rng = Random.default_rng(), pinit = rand_pinit(rs, md; rng))
    loss = select_loss(md)
    of = OptimizationFunction{true}(loss, AutoForwardDiff())
    opt_prob_params = get_opt_prob_params(md, rs)
    return OptimizationProblem(of, pinit, opt_prob_params)
end

# For an ensemble optimisation problem, creates `prob_func` function. Creates the optimisation callback.
function make_prob_func_eopt(rs::ReactionSystem, md::MeasuredData, n::Int64, traces::Vector{Vector{Float64}}, opt_prob_base;
        min_loss_thres = true_loss(md), pinits = make_pinits(rs, md, n, rng), max_loss_thres = 1e10,  rng = Random.default_rng(),
        cb_prog_print = false, cb_prog_every = 500, cb_save_trace = true)
    return function(prob, i, repeat)
        # Sets the initial condition (randomly resamples on bad starts and restarted runs).
        u0 = pinits[i]
        if (repeat > 1) || (opt_prob_base.f(u0, opt_prob_base.p) > max_loss_thres)
            u0 = rand_pinit(rs, md; rng)
        end

        # Defines the callback function.
        callback = function(opt_state, loss_val)
            # Potentially prints out progress along the optimization.
            if cb_prog_print && (opt_state.iter % cb_prog_every == 0)
                @info "Run: $i, Step: $(opt_state.iter), Loss: $(round(loss_val; sigdigits = 6)) (ADAM! ADAM! ADAM!)"
            end

            # Saves the loss function trace.
            cb_save_trace && push!(traces[i], loss_val)

            # Returns true if the loss is below the designated threshold.
            return (loss_val < min_loss_thres) || (loss_val > max_loss_thres)
        end

        # Returns the updated problem.
        return remake(prob; u0, callback)
    end
end

# For an ensemble optimisation problem, creates `output_func` function. Handles the finishing LBFGS runs.
function make_output_func_eopt(traces::Vector{Vector{Float64}}, opt_prob_base, md::MeasuredData;
        max_loss_thres = 1e10, lbfgs_maxiters = 5000, rerun_print = true, final_prints = true,
        cb_save_trace = true, success_required = false)
    return function(opt_sol1, i)
        # Check that loss value passes the maximum threshold, if not, reruns the run.
        if (opt_sol1.objective > max_loss_thres)
            rerun_print && println("Run: $i, initiated rerun due reaching too high a loss function value ($(opt_sol1.objective)).")
            return (opt_sol1, true)
        end

        # Performs a fine-tuning run using LBFGS.
        callback = function(opt_state, loss_val)
            cb_save_trace && push!(traces[i], loss_val)
            return false
        end
        opt_prob_finish = remake(opt_prob_base; u0 = opt_sol1.u, callback)
        opt_sol2 = solve(opt_prob_finish, OptimizationOptimJL.LBFGS(linesearch = LineSearches.BackTracking()); maxiters = lbfgs_maxiters)

        # Bundles the problem function and the trajectory together.
        eos = ExpandedOptimizationSolution(opt_sol1, opt_sol2, traces[i]; final_prints)

        # If we requires the final loss to be better than the true loss (should always be possible),
        # initiates a re-run.
        if success_required && objective(eos) > true_loss(md)
            rerun_print && println("Run did not surpass the true loss value ($(objective(eos)) > $(true_loss(md))).")
            return (eos, true)
        end

        return eos, false
    end
end

# For a given model and measured data, generates the n initial parameter starting values.
function make_pinits(rs::ReactionSystem, md::MeasuredData, n::Int64, rng)
    lb, ub = get_bounds(rs, md)
    samples = QuasiMonteCarlo.sample(n, lb, ub, LatinHypercubeSample(rng))
    return [[samples[:,idx]; rand_nn_p_vals(rs; rng)] for idx in 1:n]
end

# Randomly generates a single, valid, initial parameter set guess.
function rand_pinit(rs::ReactionSystem, md::MeasuredData; rng = Random.default_rng())
    lb, ub = get_bounds(rs, md)
    regulars_p_init = [rand(rng, Uniform(l, u)) for (l, u) in zip(lb, ub)]
    nn_p_init = rand_nn_p_vals(rs; rng)
    return vcat(regulars_p_init, nn_p_init)
end


### Loss Function Definition ###

# For the special case of only a single condition (frequent), uses a special loss function
# (which uses less complicated architecture, in case this matters for e.g. AD).
select_loss(md) = (numconditions(md) == 1) ? loss_single : loss_ensemble

# The loss function for a parameter set when the data contains only a single condition.
function loss_single(p, (oprob_base, measurement_ts, measurement_us, measured_sps, set_ps, u_getters, bounds, oprobs_base, set_pss, monotonicity, ps_weight_template, σ))
    p_out_of_bounds(p, bounds) && return Inf
    p = enforce_nn_monotonicity(p, monotonicity, ps_weight_template)
    p = set_ps(oprob_base, p)
    new_oprob = remake(oprob_base; p)
    new_osol = solve(new_oprob; saveat = measurement_ts, verbose = false, maxiters = 10000)
    return traj_loss(new_osol, measured_sps, measurement_us[1], u_getters, σ)
end

# The loss function for a parameter set when the data contains multiple conditions.
function loss_ensemble(p, (oprob_base, measurement_ts, measurement_us, measured_sps, set_ps, u_getters, bounds, oprobs_base, set_pss, monotonicity, ps_weight_template, σ))
    p_out_of_bounds(p, bounds) && return Inf
    p = enforce_nn_monotonicity(p, monotonicity, ps_weight_template)
    new_oprobs = [remake(oprob; p = psetter(oprob, p)) for (oprob, psetter) in zip(oprobs_base, set_pss)]
    prob_func = make_prob_func_eloss(new_oprobs, measured_sps)
    output_func = make_output_func_eloss(measurement_us, measured_sps, u_getters, σ)
    new_eprob = EnsembleProblem(new_oprobs[1]; prob_func, output_func, safetycopy = false)
    ensemble_loss = solve(new_eprob; trajectories = length(measurement_us), saveat = measurement_ts,
        verbose = false, maxiters = 10000)
    return sum(ensemble_loss)
end

# Checks if a parameter set is out of bounds.
p_out_of_bounds(p, bounds) = !all(lb <= pval <= up for (lb,pval,up) in zip(bounds[1], p, bounds[2]))

# Updates the optimisation parameter vector to ensure monotonicity for fitted neural networks.
# In practise, applied the softplus function to weight parameters.
function enforce_nn_monotonicity(p, monotonicity, ps_weight_template)
    if monotonicity == :increasing
        return [update ? Lux.softplus(p[idx]) : p[idx] for (idx,update) in enumerate(ps_weight_template)]
    elseif monotonicity == :decreasing
        return [update ? -Lux.softplus(p[idx]) : p[idx] for (idx,update) in enumerate(ps_weight_template)]
    end
    return p
end

# For a simulated trajectory, computes its loss for a given data set. Here, `measured_us`,
# is the measurement for a single condition.
# Computes the likelihood of each measurement. Assumes normally distributed noise with std = μ*σ.
# Because Optimization finds minimum, we here return the negative log-likelihood.
function traj_loss(new_osol, measured_sps, measured_us, u_getters, σ)
    SciMLBase.successful_retcode(new_osol) || return Inf
    return sum(
        sum(log(sqrt(2π)*σ) + log(max(1e-12,x_sim)) + 0.5abs2((max(1e-12,x_sim) - x_mes)/(σ*max(1e-12,x_sim)))
            for (x_sim,x_mes) in zip(u_getters[sp](new_osol), measured_us[sp]))
        for sp in measured_sps
    )
end

# When fitting a model to data with multiple condition, and ensemble simulation across all
# conditions is performed. This function returns the function initialising the individual problems.
function make_prob_func_eloss(new_oprobs, measured_sps)
    return function(prob, i, repeat)
        return new_oprobs[i]
    end
end

# When fitting a model to data with multiple condition, and ensemble simulation across all
# conditions is performed. This function evaluates a single trajectory and returns its loss.
# The first output argument is the loss, the second indicates that we will not resimulate the problem.
function make_output_func_eloss(measurement_us, measured_sps, u_getters, σ)
    return (sol, cond) -> (traj_loss(sol, measured_sps, measurement_us[cond], u_getters, σ), false)
end

### Optimisation Problem Solving ###

# Fits a model to data.
function fit_model(opt_prob, n; optalg = Adam(5e-3), ealg = EnsembleThreads(), maxiters = 10000)
    return solve(opt_prob, optalg, ealg; trajectories = n, maxiters)
end

### Optimization Problem/Solution Basic Handling ###

# For an (ensemble or normal) OptimizationSolution, find the objective value.
objective(opt_sol::ExpandedOptimizationSolution) = opt_sol.objective
objective(opt_sol::SciMLBase.OptimizationSolution) = opt_sol.objective
objective(opt_sols::SciMLBase.Union{SciMLBase.EnsembleSolution,Vector{ExpandedOptimizationSolution}}) =
    minimum(opt_sol.objective for opt_sol in opt_sols)

# For an (ensemble or normal) OptimizationSolution, find the best solution.
best_u(opt_sol::ExpandedOptimizationSolution) = opt_sol.u
best_u(opt_sol::SciMLBase.OptimizationSolution) = opt_sol.u
best_u(opt_sols::Union{SciMLBase.EnsembleSolution,Vector{ExpandedOptimizationSolution}}) =
    opt_sols[argmin([opt_sol.objective for opt_sol in opt_sols])].u

# For an (ensemble or normal) ExpandedOptimizationSolution, find the best solution when tuning is discounted.
best_u1(opt_sol::ExpandedOptimizationSolution) = opt_sol.opt_sol1.u
best_u1(opt_sols::Union{SciMLBase.EnsembleSolution,Vector{ExpandedOptimizationSolution}}) =
    opt_sols[argmin([opt_sol.opt_sol1.objective for opt_sol in opt_sols])].opt_sol1.u

# Gets The optimization parameters associated with an (ensemble/) `OptimizationProblem`.
get_opt_params(opt_sol::ExpandedOptimizationSolution) = MTK.SymbolicIndexingInterface.parameter_values(opt_sol.opt_sol1)
get_opt_params(opt_sol::SciMLBase.OptimizationSolution) = MTK.SymbolicIndexingInterface.parameter_values(opt_sol)
get_opt_params(opt_sols::SciMLBase.EnsembleSolution) = get_opt_params(opt_sols[1])
get_opt_params(opt_sols::Vector{ExpandedOptimizationSolution}) = get_opt_params(opt_sols[1])

# Gets the optimization problem from an `OptimizationProblem`/ensemble of `OptimizationProblem`s.
get_optprob(opt_prob::OptimizationProblem) = opt_prob
get_optprob(opt_probs::Optimization.EnsembleProblem) = opt_probs.prob

### Optimisation Solution Evaluations ###

# Summarises the results of the parameter fitting.
# Work in progress, can be improved to provide more information.
function summarise_optimisation(opt_sol, md::MeasuredData; print_obj = false)
    # Prints successful run information.
    successful_runs = count(obj < true_loss(md) for obj in [sol.objective for sol in opt_sol])
    if print_obj
        println("Loss thres: $(round(true_loss(md); sigdigits = 4))")
        loss_vals = [sol.objective for sol in opt_sol]
        foreach(lv -> println("Loss:       $(round(lv; sigdigits = 4))"), sort(loss_vals))
    end
    println("Successful optimisation runs: $(successful_runs)/$(length(opt_sol)) ($(round(100*successful_runs/length(opt_sol); sigdigits = 3))%)")

    # Checks the efficiency of the LBFGS follow-up optimisation.
    lbfgs_improvements = count(sol.objective < sol.opt_sol1.objective for sol in opt_sol)
    println("Optimisation runs improved by LBFGS: $(lbfgs_improvements)/$(length(opt_sol)) ($(round(100*lbfgs_improvements/length(opt_sol); sigdigits = 3))%)")
end

## Utility ###

# For a fitted parameter set, returns the `ODEProblem` corresponding to the fitted parameters.
function get_fitted_oprob(opt_sol; cond = 1)
    _, _, _, _, set_ps, _, _, oprob_bases, _ = get_opt_params(opt_sol)
    p = set_ps(oprob_bases[cond], best_u(opt_sol))
    return remake(oprob_bases[cond]; p)
end

# For a fitted parameter set, simulates the model for that parameter set and returns the simulation.
get_fitted_osol(opt_sol; cond = 1, kwargs...) = solve(get_fitted_oprob(opt_sol; cond), kwargs...)

# Retrieves the parameters that are fitted as part of the optimisation process.
function fitted_parameters(rs)
    psyms_reg = reg_parameters(rs)
    psyms_nn = vcat([MTK.wrap(p) for p in nn_vals_parameters(rs)]...)
    return [psyms_reg; psyms_nn]
end
