# A modified SIR model where the infection rate depends on both S and I.
# In the paper we absorb I into the rate, here we leave it as a reactant (but also in the
# rate for modulation).

# Fetches the parameter ranges and the models for the modified SIR model.
function get_modified_sir(ps_true, rng;
        chain = Lux.Chain(
            Lux.Dense(2 => 3, Lux.softplus, use_bias = false),
            Lux.Dense(3 => 3, Lux.softplus, use_bias = false),
            Lux.Dense(3 => 1, Lux.softplus, use_bias = false)
        )
    )
    nn_symb, nn = make_nn(chain, rng)
    return (
        get_modified_sir_ranges(),
        (
            make_modified_sir_known(ps_true),
            make_modified_sir_params(),
            make_modified_sir_nn(nn, nn_symb),
            make_modified_sir_wrong()
        )
    )
end

# Fetches the parameter ranges for the modified sir model.
get_modified_sir_ranges() = [
    :β => (-40.0, 10.0),
    :γ => (-40.0, 10.0),
    :α₁ => (-20.0, 20.0),
    :α₂ => (-20.0, 20.0)
]

# The model where the rate function is fully known (the only unknown parameter is γ).
function make_modified_sir_known(params)
    # Creates the model.
    params = make_dict(params)
    β = params[:β]
    γ = params[:γ]
    α₁ = params[:α₁]
    α₂ = params[:α₂]
    rs = @reaction_network known begin
        (10^($β)) * (S^($α₁)) * (I^($α₂)), S + I --> 2I
        10^γ, I --> R
    end

    # Sets the metadata.
    MTK.get_metadata(rs)[:rates_of_interest] = [reactions(rs)[1].rate]
    MTK.get_metadata(rs)[:dependent_sps] = [[:S, :I]]

    return rs
end

# The model where the rate function is parametrised (with β, α₁, and α₂). Only other unknown
# parameter is γ.
function make_modified_sir_params()
    # Creates the model.
    rs = @reaction_network params begin
        (10^β) * (S^α₁) * (I^α₂), S + I --> 2I
        10^γ, I --> R
    end

    # Sets the metadata.
    MTK.get_metadata(rs)[:rates_of_interest] = [reactions(rs)[1].rate]
    MTK.get_metadata(rs)[:dependent_sps] = [[:S, :I]]

    return rs
end

# The model where the rate is fitted as a neural network. Only other unknown parameter is γ.
function make_modified_sir_nn(nn, nn_symb)
    # Creates the model.
    rs = @reaction_network nn begin
        $nn(S,I), S + I --> 2I
        10^γ, I --> R
    end

    # Sets the metadata.
    MTK.get_metadata(rs)[:rates_of_interest] = [reactions(rs)[1].rate]
    MTK.get_metadata(rs)[:nns] = [nn_symb]
    MTK.get_metadata(rs)[:dependent_sps] = [[:S, :I]]

    return rs
end

# The model where the rate function is a simple parameter β (yielding the wrong model). The
# only unknown parameter is γ.
function make_modified_sir_wrong()
    # Creates the model.
    rs = @reaction_network wrong begin
        10^β, S + I --> 2I
        10^γ, I --> R
    end

    # Sets the metadata.
    MTK.get_metadata(rs)[:rates_of_interest] = [reactions(rs)[1].rate]
    MTK.get_metadata(rs)[:dependent_sps] = [[]]

    return rs
end
