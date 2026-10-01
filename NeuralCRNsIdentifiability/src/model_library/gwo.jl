# A version of the goodwind oscillator model. Some parameters are reduced to achieve structural identifiability.

# Fetches the parameter ranges and the models for the simple self activation loop model.
function get_gwo(ps_true, rng;
        chain = Lux.Chain(
            Lux.Dense(1 => 2, Lux.softplus, use_bias = false),
            Lux.Dense(2 => 2, Lux.softplus, use_bias = false),
            Lux.Dense(2 => 1, Lux.softplus, use_bias = false)
        )
    )
    nn_symb, nn = make_nn(chain, rng)
    return (
        get_gwo_ranges(),
        (
            make_gwo_known(ps_true),
            make_gwo_params(ps_true),
            make_gwo_nn(ps_true, nn, nn_symb)
        )
    )
end

# Fetches the parameter ranges and only the UDE version of the model.
function get_gwo_nn(ps_true, rng;
        chain = Lux.Chain(
            Lux.Dense(1 => 2, Lux.softplus, use_bias = false),
            Lux.Dense(2 => 2, Lux.softplus, use_bias = false),
            Lux.Dense(2 => 1, Lux.softplus, use_bias = false)
        )
    )
    nn_symb, nn = make_nn(chain, rng)
    return (
        get_gwo_ranges(),
        (
            make_gwo_nn(ps_true, nn, nn_symb),
        )
    )
end

# Fetches the parameter ranges for the goodwind oscillator model.
get_gwo_ranges() = [
    :v => (-2.0, 1.0),
    :K => (-2.0, 1.0),
    :n => (1.0, 5.0),
    :pₑ => (-2.0, 1.0),
    :pₚ => (-2.0, 1.0),
    :dₑ => (-2.0, 1.0),
    :dₘ => (-2.0, 1.0),
    :dₚ => (-2.0, 1.0)
]

# The model where the rate function is fully known (the only unknown parameter is dₚ).
function make_gwo_known(params)
    # Creates the model.
    params = make_dict(params)
    v = params[:v]
    K = params[:K]
    n = params[:n]
    pₑ = params[:pₑ]
    pₚ = params[:pₚ]
    dₑ = params[:dₑ]
    dₘ = params[:dₘ]
    rs = @reaction_network known begin
        (hillr(P, 10^($v), 10^($K), $n), 10^($dₘ)), 0 <--> M
        (10^($pₑ)*M,10^($dₑ)), 0 <--> E
        (10^($pₚ)*E,10^dₚ), 0 <--> P
    end

    # Sets the metadata.
    MTK.get_metadata(rs)[:rates_of_interest] = [reactions(rs)[1].rate]
    MTK.get_metadata(rs)[:dependent_sps] = [[:P]]

    return rs
end

# The model where the rate function is parametrised. Only other unknown parameter is dₚ.
function make_gwo_params(params)
    # Creates the model.
    params = make_dict(params)
    pₑ = params[:pₑ]
    pₚ = params[:pₚ]
    dₑ = params[:dₑ]
    dₘ = params[:dₘ]
    rs = @reaction_network params begin
        (hillr(P, 10^v, 10^K, n), 10^($dₘ)), 0 <--> M
        (10^($pₑ)*M,10^($dₑ)), 0 <--> E
        (10^($pₚ)*E,10^dₚ), 0 <--> P
    end

    # Sets the metadata.
    MTK.get_metadata(rs)[:rates_of_interest] = [reactions(rs)[1].rate]
    MTK.get_metadata(rs)[:dependent_sps] = [[:P]]

    return rs
end

# The model where the rate is fitted as a neural network. Only other unknown parameter is dₚ.
function make_gwo_nn(params, nn, nn_symb)
    # Creates the model.
    params = make_dict(params)
    pₑ = params[:pₑ]
    pₚ = params[:pₚ]
    dₑ = params[:dₑ]
    dₘ = params[:dₘ]
    rs = @reaction_network nn begin
        ($nn(P), 10^($dₘ)), 0 <--> M
        (10^($pₑ)*M,10^($dₑ)), 0 <--> E
        (10^($pₚ)*E,10^dₚ), 0 <--> P
    end

    # Sets the metadata.
    MTK.get_metadata(rs)[:rates_of_interest] = [reactions(rs)[1].rate]
    MTK.get_metadata(rs)[:nns] = [nn_symb]
    MTK.get_metadata(rs)[:dependent_sps] = [[:P]]

    return rs
end



### Neural Network Constraints Models ###


# Fetches the parameter ranges and the models for the gwo model. Includes
# models where we have put constraints on the neural network.
function get_gwo_constraints(ps_true, rng;
        chain1 = Lux.Chain(
            Lux.Dense(1 => 4, Lux.softplus, use_bias = true),
            Lux.Dense(4 => 4, Lux.softplus, use_bias = true),
            Lux.Dense(4 => 1, Lux.sigmoid_fast, use_bias = true)
        ),
        chain2 = Lux.Chain(
            Lux.Dense(1 => 4, Lux.softplus, use_bias = true),
            Lux.Dense(4 => 4, Lux.softplus, use_bias = true),
            Lux.Dense(4 => 1, Lux.softplus, use_bias = true)
        ),
        chain3 = Lux.Chain(
            Lux.Dense(1 => 4, Lux.softplus, use_bias = false),
            Lux.Dense(4 => 4, Lux.softplus, use_bias = false),
            Lux.Dense(4 => 1, Lux.softplus, use_bias = false)
        )
    )
    nn_symb1, nn1 = make_nn(chain1, rng)
    nn_symb2, nn2 = make_nn(chain2, rng)
    nn_symb3, nn3 = make_nn(chain3, rng)
    return (
        get_gwo_ranges(),
        (
            make_gwo_known(ps_true),
            make_gwo_params(ps_true),
            make_gwo_nn_mondec_bounded(ps_true, nn1, nn_symb1),
            make_gwo_nn_mondec(ps_true, nn2, nn_symb2),
            make_gwo_nn(ps_true, nn3, nn_symb3)
        )
    )
end

# The model where the rates are fitted as a neural network. Only other unknown parameter is d.
# We have constrained the neural network to be monotonously decreasing and bounded.
function make_gwo_nn_mondec_bounded(params, nn, nn_symb)
    # Creates the model.
    params = make_dict(params)
    v = params[:v]
    K = params[:K]
    n = params[:n]
    pₑ = params[:pₑ]
    pₚ = params[:pₚ]
    dₑ = params[:dₑ]
    dₘ = params[:dₘ]

    # Creates the model.
    rs = @reaction_network nn_mondec_bounded begin
        ((10^($v)) * $nn(P), 10^($dₘ)), 0 <--> M
        (10^($pₑ)*M,10^($dₑ)), 0 <--> E
        (10^($pₚ)*E,10^dₚ), 0 <--> P
    end

    # Sets the metadata.
    MTK.get_metadata(rs)[:rates_of_interest] = [reactions(rs)[1].rate]
    MTK.get_metadata(rs)[:nns] = [nn_symb]
    MTK.get_metadata(rs)[:nn_transforms] = [y -> (10^(v)) * y]
    MTK.get_metadata(rs)[:dependent_sps] = [[:P]]

    return rs
end

# The model where the rates are fitted as a neural network. Only other unknown parameter is d.
# We have constrained the neural network to be monotonously decreasing.
function make_gwo_nn_mondec(params, nn, nn_symb)
    # Creates the model.
    params = make_dict(params)
    v = params[:v]
    K = params[:K]
    n = params[:n]
    pₑ = params[:pₑ]
    pₚ = params[:pₚ]
    dₑ = params[:dₑ]
    dₘ = params[:dₘ]
    rs = @reaction_network nn_mondec begin
        ($nn(P), 10^($dₘ)), 0 <--> M
        (10^($pₑ)*M,10^($dₑ)), 0 <--> E
        (10^($pₚ)*E,10^dₚ), 0 <--> P
    end

    # Sets the metadata.
    MTK.get_metadata(rs)[:rates_of_interest] = [reactions(rs)[1].rate]
    MTK.get_metadata(rs)[:nns] = [nn_symb]
    MTK.get_metadata(rs)[:dependent_sps] = [[:P]]

    return rs
end
