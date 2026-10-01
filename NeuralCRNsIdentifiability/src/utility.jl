# Various utility functions.

### CRN-UDE Related ###

# Checks which type of parameter a parameter is.
is_reg_parameter(p) = !is_nn_vals_parameter(p) && !is_nn_parameter(p)
function is_nn_vals_parameter(p)
    MTK.hasdefault(p) || (return false)
    symstr = String(MTK.getname(p))
    (length(symstr) >= 4) && (symstr[1:4] == "NN_p") && (return true)
    return false
end
function is_nn_parameter(p)
    MTK.hasdefault(p) || (return false)
    symstr = String(MTK.getname(p))
    (length(symstr) >= 7) && (symstr[1:7] == "nn_name") && (return true)
    return false
end

# Getters for the various types of parameters.
reg_parameters(rs::MTK.AbstractSystem) = filter(is_reg_parameter, parameters(rs))
nn_vals_parameters(rs::MTK.AbstractSystem) = filter(is_nn_vals_parameter, parameters(rs))
nn_parameters(rs::MTK.AbstractSystem) = filter(is_nn_parameter, parameters(rs))

# Checks if a `ReactionSystem` contains a neural network.
is_nn_model(rs::MTK.AbstractSystem) = any(!is_reg_parameter, parameters(rs))


### Neural Network Creation ###

# From a chain and a random number generator, generates a corresponding neural network function.
# This neural network function can then be inserted into CRN models.
# Also returns the symbolic neural network representation (used for e.g. evaluating fitted neural network values).
# `nn_idx` is used in the case of multiple neural networks (updating their parameter names to differentiate them).
function make_nn(chain, rng; nn_idx = nothing)
    # Creates the symbolic neural network.
    n_input = chain.layers[1].in_dims
    n_output = chain.layers[end].out_dims
    (n_output != 1) && error("Currently only supports 1 neural network output.")
    nn_name = make_nn_name(nn_idx)
    nn_p_name = make_p_nn_name(nn_idx)
    NN_symb, NN_ps = SymbolicNeuralNetwork(; chain, n_input , n_output, rng, nn_name, nn_p_name)

    return NN_symb, make_nn_func(n_input, NN_symb, NN_ps)
end

# Generates names for neural network symbolic parameters.
make_nn_name(nn_idx::Nothing) = :NN
make_nn_name(nn_idx::Integer) = Symbol("NN_$nn_idx")
make_p_nn_name(nn_idx::Nothing) = :NN_ps
make_p_nn_name(nn_idx::Integer) = Symbol("NN_ps_$nn_idx")

# Returns the appropriate output neural network function.
function make_nn_func(n_input, NN_symb, NN_ps)
    (n_input == 1) && return x -> NN_symb([x], NN_ps)[1]
    (n_input == 2) && return (x,y) -> NN_symb([x,y], NN_ps)[1]
    (n_input == 3) && return (x,y,z) -> NN_symb([x,y,z], NN_ps)[1]
end


### Neural Network Related ###

# Evaluates a neural network across a range of values.
function eval_nn(vals, rs::ReactionSystem, nn_idx, rs_p_vals, md::MeasuredData,
        oprob_base = nothing, set_ps = nothing)
    # Retrieves the neural network parameter values.
    if isnothing(oprob_base) || isnothing(set_ps)
        oprob_base, _, _, _, set_ps, _, _, _, _, _, _ = get_opt_prob_params(md, rs)
    end
    p = set_ps(oprob_base, rs_p_vals)
    nn_p_vals = remake(oprob_base; p).ps[nn_vals_parameters(rs)[nn_idx]]

    # Creates he neural network function, and applies it to the input values.
    nn_func_pre1 = MTK.getdefault(get_nns(rs)[nn_idx])
    nn_func_pre2(val) = nn_func_pre1(make_vec(val), nn_p_vals)[1]
    if haskey(MTK.get_metadata(rs), :nn_transforms)
        nn_func(val) = MTK.get_metadata(rs)[:nn_transforms][nn_idx](nn_func_pre2(val))
    else
        nn_func = nn_func_pre2
    end
    return nn_func.(vals)
end
make_vec(val::Float64) = [val]
make_vec(val::Vector) = val
make_vec(val::Tuple) = collect(val)

# A function which, for a symbolic neural network, generates a random set of values for its parameters.
# The `rng` and `idx` arguments` makes it reproducible and usable for ensemble optimisation initial guesses.
function rand_nn_p_vals(NN; idx = 1, rng = Random.default_rng())
    rng_idx = StableRNG(idx + rand(rng, 1:1000))
    return Vector(ComponentArray(Lux.initialparameters(rng_idx, NN.metadata[Symbolics.VariableDefaultValue].lux_model)))
end

# For a neural network, returns the full (single) vector with all neural network parameter.
# The index can be changed to generate different initial values guesses for different ensemble
# simulation runes. For a model without neural networks, returns an empty vector.
function rand_nn_p_vals(rs::ReactionSystem; idx = 1, rng = Random.default_rng())
    haskey(MTK.get_metadata(rs), :nns) || return Float64[]
    return vcat([rand_nn_p_vals(nn; idx, rng) for nn in MTK.get_metadata(rs)[:nns]]...)
end

### Parameter transformations ###

# Transformations to and from log-space.
trans_log(p::Number) = log10(p)
trans_log(ps::Vector{T}) where {T <: Number} = [trans_log(p) for p in ps]
trans_log(ps::Vector{Pair{S,T}}) where {S,T <: Number} = [p[1] => trans_log(p[2]) for p in ps]
trans_log(ps::Dict{S,T}) where {S,T <: Number} = Dict(p[1] => trans_log(p[2]) + 1 for p in ps)
trans_log_inv(p::Number) = 10^p
trans_log_inv(ps::Vector{T}) where {T <: Number} = [trans_log_inv(p) for p in ps]
trans_log_inv(ps::Vector{Pair{S,T}}) where {S,T <: Number} = [p[1] => trans_log_inv(p[2]) for p in ps]
trans_log_inv(ps::Dict{S,T}) where {S,T <: Number} = Dict(p[1] => trans_log_inv(p[2]) + 1 for p in ps)

# Transformations to and from unbounded log-space (tanh-transformation).


### Unsorted ###

# Checks if a crn model contains a neural network.
is_nn_model(rs) = !all(is_reg_parameter, parameters(rs))

# Retrieves various information stored in the `ReactionSystem` metadata.
get_rates_of_interest(rs::ReactionSystem) = rs.metadata[:rates_of_interest]
get_dependent_sps(rs::ReactionSystem) = rs.metadata[:dependent_sps]
get_nns(rs::ReactionSystem) = rs.metadata[:nns]
num_nns(rs::ReactionSystem) = length(get_nns(rs))

# Gets the length of a parameter range.
range_length(range::Tuple) = range[2] - range[1]

# Gets a `ReactionSystem`'s name (exported).`
getname(rs::ReactionSystem) = MTK.getname(rs)

# Checks if a model is "correct" (i.e. contain the true model).
is_valid(rs::ReactionSystem) = !contains(string(MTK.nameof(rs)), "wrong")

# Gets whether a certain f_index is a nn or a func.
function f_idx_is_nn(rs::ReactionSystem, f_idx::Int)
    if haskey(MTK.get_metadata(rs), :is_nn_rate)
        return MTK.get_metadata(rs)[:is_nn_rate][f_idx]
    else
        return is_nn_model(rs)
    end
end

# Converts a certain f index to a nn index.
function f_idx_to_nn_idx(rs::ReactionSystem, f_idx::Int)
    if haskey(MTK.get_metadata(rs), :is_nn_rate)
        return count(MTK.get_metadata(rs)[:is_nn_rate][1:f_idx])
    else
        return f_idx
    end
end


### Generic Stuff ###

# A function for merging keyword arguments.
function merge_kwargs(base_args...; kwargs...)
    # make NamedTuple from base_args
    base_nt = if length(base_args) == 0
        (;)
    elseif length(base_args) == 1 && isa(base_args[1], NamedTuple)
        base_args[1]
    elseif length(base_args) == 1 && isa(base_args[1], AbstractDict)
        (; collect(base_args[1])...)
    else
        (; base_args...)
    end

    # incoming kwargs as NamedTuple
    incoming_nt = (; kwargs...)

    # base should overwrite incoming -> make incoming first, then base
    return merge(incoming_nt, base_nt)
end

"""
    make_fig(plt, filename; filefolder = "Figures", filetypes = ["png,", "svg"])

Helper function which takes a figure and saves it in multiple file formats. Also returns the
figure itself (so that it is automatically displayed after `begin ... end` blocks).

Arguments:
- `plt`: The plot object to be saved.
- `filename`: The name of the file to save the figure as (without extension).
- `filefolder = "Figures"`: The folder to save the figure in. Default is "Figures".
- `filetypes = ["png", "svg"]`: The file types to save the figure as. Default is PNG and SVG.
"""
function make_fig(plt, filename; filefolder = "Figures", filetypes = ["png", "svg"], legend = :none)
    for filetype in filetypes
        plt_saved = if filetype == "png"
            plt
        else
            plot!(deepcopy(plt); formatter = :none, legend, bg = :transparent,
                bg_inside = :transparent)
        end
        savefig(plt_saved, joinpath(filefolder, "$filename.$filetype"))
    end
    return plt
end
