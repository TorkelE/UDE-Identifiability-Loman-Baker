# Contains functions for generating synthetic data from a model.

# For:
# - A model.
# - A list of initial conditions (one set of measurements is performed for each).
# - A time end point.
# - A true parameter set.
# - A list of measures species.
# - A measurement count.
# - A measurement noise function.
# Generates a synthetic dataset. This dataset is represented in a special structure.
struct MeasuredData
    rs::ReactionSystem
    u0s::Vector{Dict{Symbol,Float64}}
    tend::Float64
    ps_true::Dict{Symbol,Float64}
    ps_ranges::Dict{Symbol, Tuple{Float64,Float64}}
    measured_sps::Vector
    measurement_ts::AbstractRange{Float64}
    measurement_us::Vector{Dict{Symbol,Vector{Float64}}}
    oprobs_base::ODEProblem
    sols_true::Vector{ODESolution}
    σ::Float64

    # The constructor with just the field inputs.
    function MeasuredData(rs, u0s, tend, ps_true, ps_ranges, measured_sps, measurement_ts,
            measurement_us, oprob_base, sols_true, σ)
        new(rs, u0s, tend, ps_true, ps_ranges, measured_sps, measurement_ts, measurement_us,
            oprob_base, sols_true, σ)
    end
end

# The normal constructor. Performs simulations and stores results.
function MeasuredData(rs, u0s, tend, ps_true, sample_n, measured_sps, σ;
        rng = Random.default_rng(), ps_ranges = initialise_p_ranges(ps_true), noise_func = scaled_normal_error)
    # Error checks.
    any(v[2][1] >= v[2][2] for v in ps_ranges) && error("The second value in a parameter range must be larger than the first one.")

    # Create dictionaries of the relevant inputs.
    u0s = make_u0_dicts(u0s)
    ps_true = make_dict(ps_true)
    ps_ranges = make_dict(ps_ranges)

    # Creates the base `ODEProblem` and time stops. Generates the data from these.
    ps_relevant = relevant_ps(rs, ps_true)
    measurement_ts = range(0.0, tend, length = sample_n)
    oprob_base = ODEProblem(rs, u0s[1], tend, ps_relevant)
    noise_func = noise_func(σ; rng)
    measurements = [make_data(oprob_base, u0, measured_sps, measurement_ts, noise_func) for u0 in u0s]
    measurement_us = getindex.(measurements, 1)
    sols_true = getindex.(measurements, 2)

    # Creates the `MeasuredData` structure.
    MeasuredData(rs, u0s, tend, ps_true, ps_ranges, measured_sps, measurement_ts, measurement_us,
        oprob_base, sols_true, σ)
end

# A constructor which takes a "settings" dictionary (used for various HPC scans).
function MeasuredData(settings::Dict{Symbol,Any})
    rng = StableRNG(settings[:seed])
    ps_ranges, models = settings[:retrieval_function](settings[:ps_true], rng)
    MeasuredData(models[1], settings[:u0s], settings[:tend], settings[:ps_true], settings[:sample_n],
        settings[:measured_sps], settings[:σ_std]; rng, ps_ranges, noise_func = settings[:noise_func])
end

# Loads a measured data set that have been saved to file(s).
function MeasuredData(file_path::String)
    # If stored as a ".jls" file, it is a stored Julia object we can load directly.
    file_path[end:end-3] == ".jls" || return deserialize(file_path)

    # Loads data that has been written manually to file.
    error("Loading `MeasuredData` from non-Julia serialisation files is currently not supported.")
end


### NOTE: MAKE THIS WORK WITH INPUT RNGs!!! ###
# Function for creating the measurements for an individual condition.
function make_data(oprob_base, u0, measured_sps, measurement_ts, noise_func)
    # Simulates the model.
    oprob = remake(oprob_base; u0)
    sol_true_samples = solve(oprob, Rodas5P(); saveat = measurement_ts, abstol = 1e-8, reltol = 1e-8)
    sol_true = solve(oprob, Rodas5P(); abstol = 1e-8, reltol = 1e-8)

    # Creates the measurement output. Goes through it an applied the noise function to each value.
    # Returns these as well as the true simulation.
    measurements = Dict([sps => sol_true_samples[sps] for sps in measured_sps])
    for sps in measured_sps
        for (idx,u) in enumerate(measurements[sps])
            measurements[sps][idx] = noise_func(u)
        end
    end
    return measurements, sol_true
end

# Initiates a set of parameter from their initial values.
function initialise_p_ranges(ps_true::Vector{Pair{S,T}}; factor = 50.0) where {S,T}
    Dict([p => (p_val/factor, p_val*factor) for (p, p_val) in ps_true])
end

### Noised Functions ###

# Generates a function generating a normal distribution with a certain std.
normal_error(σ; rng = Random.default_rng()) = (x -> x + rand(rng, Normal(0.0, σ)))

# Generates a function generating a normal distribution where the std scale with the mean.
scaled_normal_error(σ; rng = Random.default_rng()) = (x -> x + rand(rng, Normal(0.0, x*σ)))

### Data Subsampling ###

# Subsamples a `MeasuredData` structure, selecting less samples, conditions, and measured species.
# Crucially, the core data is preserved, so a subsampled version should always have less information.
## NOT IMPLEMENTED YET - MIGHT NOT BE NEEDED ###

### MeasuredData Accessors ###

# Counts the number of measured species.
nummeasuredsps(md::MeasuredData) = length(md.measured_sps)

# Counts the number of conditions.
numconditions(md::MeasuredData) = length(md.u0s)

# Retrieves the non-measured species.
get_non_measured(md::MeasuredData) = filter(sp -> !(sp in md.measured_sps), MTK.getname.(species(md.rs)))

# Get the number of measured and non-measured species, respectively.
num_measured(md::MeasuredData) = length(md.measured_sps)
num_nonmeasured(md::MeasuredData) = length(get_non_measured(md))

### Expanded MeasuredData Accessors ###

# For a model, extracts the relevant parameters and returns their bonds as separate `lb` & `ub` vectors.
function get_bounds(rs::ReactionSystem, md::MeasuredData)
    lb = [md.ps_ranges[MTK.getname(p)][1] for p in reg_parameters(rs)]
    ub = [md.ps_ranges[MTK.getname(p)][2] for p in reg_parameters(rs)]
    return lb, ub
end

# Compute the true loss value for a `MeasuredData` structure.
function true_loss(md::MeasuredData)
    opt_prob_params = get_opt_prob_params(md, md.rs)
    ps_true = [md.ps_true[MTK.getname(p)] for p in reg_parameters(md.rs)]
    return (select_loss(md))(ps_true, opt_prob_params)
end


### Real Data ###

#  Creates a `MeasuredData` structure for a real data set.
function import_data(filename, rs, u0s, ps_true, measured_sps, σ; ps_ranges = initialise_p_ranges(ps_true))
    # Create dictionaries of the relevant inputs.
    u0s = make_u0_dicts(u0s)
    ps_true = make_dict(ps_true)
    ps_ranges = make_dict(ps_ranges)

    # Reads data.
    data = CSV.File(filename)
    # measurement_ts = Float64.(data.columns[1].column)
    # measurement_ts = measurement_ts[1]:1.0:measurement_ts[end]
    measurement_ts = 1.0:length(data)
    tend = measurement_ts[end]
    measurement_us = [Dict(measured_sps[1] => Float64.(data.columns[2].column))]

    # Finalises the `MeasuredData` structure.
    ps_relevant = relevant_ps(rs, ps_true)
    oprob_base = ODEProblem(rs, u0s[1], tend, ps_relevant)
    sols_true = [solve(oprob_base, Rodas5P(); abstol = 1e-8, reltol = 1e-8)]

    MeasuredData(rs, u0s, tend, ps_true, ps_ranges, measured_sps, measurement_ts, measurement_us,
        oprob_base, sols_true, σ)
end

### Unsorted Functions ###

# Creates the parameters to the optimization problem (not the ones we are actually optimising).
function get_opt_prob_params(md::MeasuredData, rs::ReactionSystem; tend = md.tend)
    psyms_optim = fitted_parameters(rs)
    ps_base = relevant_ps(rs, md.ps_true)
    oprob_base = ODEProblem(rs, md.u0s[1], tend, ps_base)
    oprob_bases = [ODEProblem(rs, u0, tend, ps_base) for u0 in md.u0s]
    set_ps = MTK.setp_oop(oprob_base, psyms_optim)
    bounds = get_bounds(rs, md)
    set_pss = [MTK.setp_oop(oprob, psyms_optim) for oprob in oprob_bases] # Have to be added due to a bug where setp will update initial conditions.
    get_us = make_u_getters(oprob_base, md)
    measurement_ts = filter(t -> t <= tend, md.measurement_ts)
    measurement_us = [
        Dict([ sp => md.measurement_us[cond][sp][1:length(measurement_ts)]
            for sp in md.measured_sps])
        for cond in eachindex(md.measurement_us)
    ]
    monotonicity = :false
    occursin("moninc", String(getname(rs))) && (monotonicity = :increasing)
    occursin("mondec", String(getname(rs))) && (monotonicity = :decreasing)
    ps_weight_template = if monotonicity != :false
        NN = only(MTK.get_metadata(rs)[:nns])
        nn_p_labels = ComponentArrays.labels(ComponentArray(Lux.initialparameters(Random.default_rng(), NN.metadata[Symbolics.VariableDefaultValue].lux_model)))
        nn_ps_weight_template = [occursin("weight", label) for label in nn_p_labels]
        vcat(fill(false, length(reg_parameters(rs))), nn_ps_weight_template)
    else
        fill(false, length(ps_base))
    end
    return (oprob_base, measurement_ts, measurement_us, md.measured_sps, set_ps, get_us, bounds, oprob_bases, set_pss, monotonicity, ps_weight_template, md.σ)
end

# From a base `ODEProblem` and some measured data,` generates a getu function for use
# within the optimisation loss function.
function make_u_getters(oprob::ODEProblem, md::MeasuredData)
    Dict([sp => MTK.getu(oprob, sp) for sp in md.measured_sps])
end

### Utility Functions ###

# Converts potential maps to dictionaries.
make_dict(map::Dict) = map
make_dict(map::Vector{Pair{S,T}}) where {S,T} = Dict(map)

# Converts initial conditions to vector of dictionaries.
make_u0_dicts(u0s::Dict{S,T}) where {S,T} = [u0s]
make_u0_dicts(u0s::Vector{Pair{S,T}}) where {S,T} = [make_dict(u0s)]
make_u0_dicts(u0s::Vector{Dict{S,T}}) where {S,T} = u0s
make_u0_dicts(u0s::Vector{Vector{Pair{S,T}}}) where {S,T} = make_dict.(u0s)

# Filters out parameters that are not part of the given system.
relevant_ps(rs, ps) = filter(ppair -> ppair[1] ∈ MTK.getname.(parameters(rs)), ps)
