## Common Likelihood Profiling methods ###

# A template with the optimisation profiler method settings.
# In practise, this (and other LikelihoodProfiler.jl-based approaches) were not used.
struct OptimisationPLMethod
    m::Int64
    optimizer
    maxiters::Int64

    function OptimisationPLMethod(;m = 100, optimizer = NLopt.LN_BOBYQA(), maxiters = 5e4)
        new(m, optimizer, maxiters)
    end
end

# A template with the integration profiler method settings.
# In practise, this (and other LikelihoodProfiler.jl-based approaches) were not used.
struct IntegrationPLMethod
    m::Int64
    integrator
    maxiters::Int64

    function IntegrationPLMethod(;m = 100, integrator = Tsit5(), maxiters = 5e4)
        new(m, integrator, maxiters)
    end
end

# A template with the custom profiler method settings.
struct CustomPLMethod
    m::Int64
    n_starts::Int64
    adam_maxiters::Int64
    inf_stop::Bool
    max_thres::Float64
    p_range::Union{Nothing, Tuple{Float64, Float64}}

    function CustomPLMethod(;p_range = nothing, m = 10, n_starts = 5, adam_maxiters = 10000,
            inf_stop = true, max_thres = 1e10)
        new(m, n_starts, adam_maxiters, inf_stop, max_thres, p_range)
    end
end

# A structure which stores the relevant likelihood profile information.
struct LikelihoodProfile
    pars::Vector{Vector{Float64}}
    x::Vector{Float64}
    obj::Vector{Float64}
    retcodes::Tuple{Symbol, Symbol}
    mle_idx::Int64
end


### Likelihood Profiling via LikelihoodProfiler ###

# Computes all likelihood profiles for a given models (not related to neural networks).
function compute_lprofiles(method, rs, opt_prob, opt_sol, md; maxiters = 5e4, print = false, ealg = EnsembleThreads(), kwargs...)
    return Dict(p => compute_lprofile(method, MTK.getname(p), rs, opt_prob, opt_sol, md; print, ealg, kwargs...)
        for p in reg_parameters(rs))
end

# Computes the liklihood profile for a specific parameter.
function compute_lprofile(method::Union{OptimisationPLMethod,IntegrationPLMethod}, p::Symbol, rs::ReactionSystem, opt_prob, opt_sol, md; print = false)
    # Error check.
    p_idx = findfirst(p .== MTK.getname.(fitted_parameters(rs)))
    isnothing(p_idx) && error("Likelihood parameter ($p) not found among model parameters ($parameters(rs)).")
    (numparams(rs) == 1) && error("Likelihood profiling currently not supported for models with only a single parameter.")

    # Creates the `PLProblem` and computes the profile.
    profile_ranges = [(v - 1e-4, v + 1e-4) for v in best_u(opt_sol)]
    profile_ranges[p_idx] = md.ps_ranges[p]
    plprob = PLProblem(get_optprob(opt_prob), best_u(opt_sol), profile_ranges; threshold = Inf)
    profiler = make_method(method, range_length(md.ps_ranges[p]))
    pl_sol = profile(plprob, profiler; idxs = [p_idx], maxiters = method.maxiters)[1]

    # Uses a rough method to compute the maximum likelihood estimate index.
    mle_idx = findfirst(pl_sol.obj .== minimum(pl_sol.obj))
    (minimum(pl_sol.obj) != objective(opt_sol)) && @warn "The minimum likelihood value in the profile does not match the objective value of the optimisation solution."

    # Checks the output and then returns it.
    (pl_sol.retcodes[1] == :MaxIters) && @warn "Likelihood profiling reach MaxIters (lower bound)."
    (pl_sol.retcodes[2] == :MaxIters) && @warn "Likelihood profiling reach MaxIters (upper bound)."
    return LikelihoodProfile(pl_sol.pars, pl_sol.x, pl_sol.obj, pl_sol.retcodes, mle_idx)
end

# From a profiling template, generates the corresponding profiling method.
function make_method(method::OptimisationPLMethod, range_length)
    initial_step = range_length / method.m
    stepper = FixedStep(; initial_step)
    optimizer_opts = (; xtol_abs = 1e-7,)
    return OptimizationProfiler(; optimizer = method.optimizer, optimizer_opts, stepper)
end
function make_method(method::IntegrationPLMethod, range_length)
    dtmax = range_length / method.m
    integrator_opts = (; dtmax,)
    return IntegrationProfiler(; integrator = method.integrator, integrator_opts, matrix_type = :hessian)
end

### Custom-written Likelihood Profiling ###

# Computes the liklihood profile for a specific parameter.
# Uses customised algorithm.
# The `opt_prob` argument is not actually used here.
function compute_lprofile(method::CustomPLMethod, p::Symbol, rs::ReactionSystem, opt_prob, opt_sol, md;
        rng = Random.default_rng(), print = false, dynamic_max_thres = 10.0, p_range = (md.ps_ranges[p][1], md.ps_ranges[p][2]), ealg = EnsembleDistributed())
    # Prints the profiling information.
    print && println("\nComputing likelihood profile for parameter $p.")

    # Fetches method settings. Potentially adapts for the dynamic maximum threshold.
    m, inf_stop, max_thres = (method.m, method.inf_stop, method.max_thres)
    max_thres = min(max_thres, objective(opt_sol) + dynamic_max_thres)

    # Error check.
    p_idx = findfirst(p .== MTK.getname.(fitted_parameters(rs)))
    isnothing(p_idx) && error("Likelihood parameter ($p) not found among model parameters ($parameters(rs)).")

    # Creates the grid over which we compute the profile.
    isnothing(method.p_range) || (p_range = method.p_range)
    prof_grid = collect(range(p_range[1]; stop = p_range[2], length = m))
    p_start = best_u(opt_sol)[p_idx]
    start_idx = findlast(p_start .> prof_grid)
    isnothing(start_idx) && (start_idx = 0)
    start_idx = start_idx + 1
    insert!(prof_grid, start_idx, p_start)

    # Initialises the storage of the profiling info.
    prof_vals = fill(Inf, length(prof_grid))
    prof_sols = Vector{Vector{Float64}}(undef, length(prof_grid))
    prof_vals[start_idx] = objective(opt_sol)
    prof_sols[start_idx] = deleteat!(deepcopy(best_u(opt_sol)), p_idx)

    # Solves the optimisation problem along each grid point. If the new step is better than
    # the last one, potentially reruns that step using the new step as minimum threshold.
    p_symbolic = fitted_parameters(rs)[p_idx]
    idx = start_idx + 1
    while (idx <= length(prof_grid))
        print && println("Forward direction: Computing likelihood profile for at idx $idx.")
        rerun_thres = minimum(prof_vals[idx:end])
        lprofile_step!(idx, idx - 1, idx + 1, rs, md, p_symbolic, prof_grid, prof_sols, prof_vals, rng,
            method, rerun_thres, ealg)
        if inf_stop && (prof_vals[idx] >= max_thres) && (prof_vals[idx-1] >= max_thres)
            break
        end
        idx_relative = findfirst(prof_vals[start_idx+1:end][i] - 0.2 * max(idx - start_idx - i,0) > prof_vals[idx] for i in eachindex(prof_vals[start_idx+1:end]))
        isnothing(idx_relative) && break
        if start_idx + idx_relative < idx
            println("Stepping backwards: $(prof_vals[start_idx + idx_relative]) > $(prof_vals[idx])")
        end
        idx = start_idx + idx_relative
    end
    idx = start_idx - 1
    while idx >= 1
        print && println("Backwards direction: Computing likelihood profile for at idx $idx.")
        rerun_thres = minimum(prof_vals[1:idx])
        lprofile_step!(idx, idx + 1, idx - 1, rs, md, p_symbolic, prof_grid, prof_sols, prof_vals, rng,
            method, rerun_thres, ealg)
        if inf_stop
            if (prof_vals[idx] >= max_thres) && (prof_vals[idx + 1] >= max_thres)
                break
            end
        end
        idx_new = findlast(prof_vals[1:start_idx-1] .> prof_vals[idx] + 0.2)
        if !isnothing(idx_new) && (idx_new > idx)
            println("Stepping backwards: $(prof_vals[idx_new]) > $(prof_vals[idx])")
        end

        idx = findlast(prof_vals[1:start_idx-1] .> prof_vals[idx])
        isnothing(idx) && break
    end

    # performs post-processing identifications and smoothing of bad points.
    try_count = 0
    idx_up_count = zeros(length(prof_vals))
    elap = @elapsed while try_count < 100
        stand_out_vals = filter(idx -> stand_out(prof_vals, idx), 7:length(prof_vals)-6)
        stand_out_vals = filter(idx -> idx_up_count[idx] < 6, stand_out_vals)
        println("\nPostprocessing step $try_count: Found $(length(stand_out_vals)) stand-out values.")
        println(stand_out_vals)
        isempty(stand_out_vals) && break
        println(prof_vals[stand_out_vals])
        idx = stand_out_vals[argmax(prof_vals[stand_out_vals])]
        println(idx)
        idx_up_count[idx] += 1
        if idx > start_idx
            rerun_thres = minimum(prof_vals[idx:end])
            lprofile_step!(idx, idx - 1, idx + 1, rs, md, p_symbolic, prof_grid, prof_sols, prof_vals, rng,
                method, rerun_thres, ealg)
        else
            rerun_thres = minimum(prof_vals[1:idx])
            lprofile_step!(idx, idx + 1, idx - 1, rs, md, p_symbolic, prof_grid, prof_sols, prof_vals, rng,
                method, rerun_thres, ealg)
        end
        try_count += 1
    end
    print && println("Number of postprocessing smoothing steps: $try_count. Total time: $(round(elap; sigdigits = 3)) seconds.")

    # Uses a rough method to compute the maximum likelihood estimate index.
    mle_idx = findfirst(prof_vals .== minimum(prof_vals))
    (minimum(prof_vals) != objective(opt_sol)) && @warn "The minimum likelihood value in the profile does not match the objective value of the optimisation solution."

    return LikelihoodProfile(prof_sols, prof_grid, prof_vals, (:NonIdentifiable, :NonIdentifiable), mle_idx)
end

# Performs a step of the custom likelihood profiler.
function lprofile_step!(idx_new::Int64, idx_prev::Int64, idx_next::Int64, rs::ReactionSystem, md::MeasuredData,
        p_symbolic, prof_grid, prof_sols, prof_vals, rng, method, rerun_thres, ealg)
    # Extracts the method settings.
    n_starts, maxiters, max_thres = (method.n_starts, method.adam_maxiters, method.max_thres)

    # Creates an ensemble optimisation problem.
    fixed_p = (p_symbolic => prof_grid[idx_new])
    p_seed = prof_sols[idx_prev]
    eoptprob = lprofile_make_eoptprob(rs, md, n_starts, fixed_p, p_seed; rng, rerun_thres, max_thres)
    single_optprob = OptimizationProblem(eoptprob.prob.f, p_seed, eoptprob.prob.p)

    # Potentially creates a second single optimisation problem for the next index (if it is better).
    single_optprob_2 = if (1 <= idx_next <= length(prof_grid)) && (prof_vals[idx_next] < prof_vals[idx_new])
        p_seed2 = prof_sols[idx_next]
        OptimizationProblem(eoptprob.prob.f, p_seed2, eoptprob.prob.p)
    else
        nothing
    end

    # Solves the optimisation problem (potentially trivially) and extracts the results.
    if length(eoptprob.prob.u0) > 0
        # Solves the optimization problems.
        lp_eopt_sol = solve(eoptprob, Adam(5e-3), ealg;
            maxiters, trajectories = n_starts)
        single_lp_opt_sol = solve(single_optprob, OptimizationOptimJL.LBFGS(linesearch = LineSearches.BackTracking()); maxiters = 1000)

        # If we used a second forward seed, solves that optimisation problem and utilises it.
        if !isnothing(single_optprob_2)
            single_lp_opt_sol_2 = solve(single_optprob_2, OptimizationOptimJL.LBFGS(linesearch = LineSearches.BackTracking()); maxiters = 1000)
            if objective(single_lp_opt_sol_2) < objective(single_lp_opt_sol)
                single_lp_opt_sol = single_lp_opt_sol_2
            end
        end

        # Finds the correct objective value and solution to store.
        if objective(lp_eopt_sol) < objective(single_lp_opt_sol)
            prof_vals[idx_new] = objective(lp_eopt_sol)
            prof_sols[idx_new] = best_u(lp_eopt_sol)
        else
            prof_vals[idx_new] = objective(single_lp_opt_sol)
            prof_sols[idx_new] = best_u(single_lp_opt_sol)
        end
    else
        prof_vals[idx_new] = eoptprob.prob.f([], eoptprob.prob.p)
        prof_sols[idx_new] = []
    end
    (prof_vals[idx_new] > max_thres) && (prof_vals[idx_new] = Inf)
    # markercolor = fill(1, 1, length(prof_vals))
    # markercolor[1, idx_new] = 2
    # display(plot(prof_grid, -prof_vals; markershape = :circle, ms = 2.0), markercolor, label = "", xlimit = (prof_grid[1], prof_grid[end]))
end

# Checks if a value stands out.
function stand_out(prof_vals, idx)
    (prof_vals[idx] > prof_vals[idx + 1] > prof_vals[idx + 2] > prof_vals[idx + 3]) &&
        (prof_vals[idx] > prof_vals[idx - 1] > prof_vals[idx - 2] > prof_vals[idx - 3]) && return false
    return prof_vals[idx] > (max(prof_vals[idx - 1], prof_vals[idx+1]) + 0.1)
end

### Manual LikelihoodProfiler Problem-generating Functions ###

# Makes an ensemble optimisation problem with one parameter value fixed.
function lprofile_make_eoptprob(rs::ReactionSystem, md::MeasuredData, n::Int64, fixed_p, p_seed;
        rng = Random.default_rng(), rerun_thres = Inf, max_thres = Inf)
    opt_prob_base = lprofile_make_optprob(fixed_p, p_seed, rs, md; rng)
    prob_func = lprofile_make_prob_func_eopt(p_seed, fixed_p, rs, md, opt_prob_base, n; max_thres, rng)
    output_func = lprofile_make_output_func_eopt(opt_prob_base; rerun_thres)
    return EnsembleProblem(opt_prob_base; prob_func, output_func, safetycopy = false)
end

# Creates the base optimisation problem with one parameter value fixed.
function lprofile_make_optprob(fixed_p, p_seed, rs, md; rng = Random.default_rng())
    opt_params = lprofile_get_opt_prob_params(fixed_p, rs, md)
    opt_prob = make_optprob(rs, md; rng)
    return remake(opt_prob; p = opt_params, u0 = p_seed)
end

# Creates the optimisation parameters for the base optimisation problem with one parameter value fixed.
function lprofile_get_opt_prob_params(fixed_p, rs, md)
    psyms_filtered = filter(!isequal(fixed_p[1]), fitted_parameters(rs))
    oprob_base, measurement_ts, measurement_us, measured_sps, _, u_getters, bounds, oprobs_base, _, monotonicity, ps_weight_template, σ =
        get_opt_prob_params(md, rs)
    oprob_base = remake(oprob_base; p = [fixed_p])
    oprobs_base = [remake(oprob; p = [fixed_p]) for oprob in oprobs_base]
    set_ps = MTK.setp_oop(oprob_base, psyms_filtered)
    set_pss = [MTK.setp_oop(oprob, psyms_filtered) for oprob in oprobs_base]
    p_idx = findfirst(isequal(fixed_p[1], p) for p in fitted_parameters(rs))
    bounds = (deleteat!(deepcopy(bounds[1]), p_idx), deleteat!(deepcopy(bounds[2]), p_idx))
    ps_weight_template = deleteat!(deepcopy(ps_weight_template), p_idx)
    return (oprob_base, measurement_ts, measurement_us, measured_sps, set_ps, u_getters, bounds, oprobs_base, set_pss, monotonicity, ps_weight_template, σ)
end

# Creates the ensemble problem initiation function for when one parameter is eliminated.
function lprofile_make_prob_func_eopt(p_seed, fixed_p, rs::ReactionSystem, md::MeasuredData,
        opt_prob_base, n::Int64; max_thres = Inf, rng = Random.default_rng())
    # Creates the new parameter initial guesses (with the seed and one parameter removed).
    p_idx = findfirst(isequal(fixed_p[1]), fitted_parameters(rs))
    seeded_pinits = make_pinits(rs, md, n-1, rng)
    seeded_pinits = [deleteat!(pinit, p_idx) for pinit in seeded_pinits]
    pinits = [[p_seed]; seeded_pinits]

    # Creates the problem function.
    return function(prob, i, repeat)
        # Sets the initial condition (randomly resamples on bad starts and restarted runs).
        # Maximally resamples 20 times.
        u0 = pinits[i]
        (repeat < 3) && for i = 1:20
            if opt_prob_base.f(u0, opt_prob_base.p) >= max_thres

                u0 = rand_pinit(rs, md; rng)
                deleteat!(u0, p_idx)
            else
                break
            end
        end

        # Defines a callback function.
        callback = (opt_state, loss_val) -> (loss_val > max_thres)

        # Returns the updated problem.
        return remake(prob; u0, callback)
    end
end

# Creates the ensemble result output function for when one parameter is eliminated.
function lprofile_make_output_func_eopt(opt_prob_base; rerun_thres = Inf, lbfgs_maxiters = 5000)
    return function(opt_sol1, i)
        # Performs a fine-tuning run using LBFGS.
        opt_prob_finish = remake(opt_prob_base; u0 = opt_sol1.u)
        opt_sol2 = solve(opt_prob_finish, OptimizationOptimJL.LBFGS(linesearch = LineSearches.BackTracking()); maxiters = lbfgs_maxiters)
        eos = ExpandedOptimizationSolution(opt_sol2, opt_sol1, [1.0]; final_prints = false)
        rerun = if (objective(eos) > rerun_thres + 0.3)
            if (rand() < 0.95)
                true
            else
                false
            end
        else
            false
        end
        return eos, rerun
    end
end
