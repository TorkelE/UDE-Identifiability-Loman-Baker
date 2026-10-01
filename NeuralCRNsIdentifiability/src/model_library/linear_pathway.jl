# A version of the linear_pathway model. This uses only 4 species and repressive interactions.
# We uses the same rate function form for all steps.

# Fetches the parameter ranges and the models for the linear pathway model.
function get_linear_pathway(ps_true, rng;
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
        get_linear_pathway_ranges(),
        (
            make_linear_pathway_known(ps_true),
            make_linear_pathway_params(),
            make_linear_pathway_nn_mondec_bounded(ps_true, nn1, nn_symb1),
            make_linear_pathway_nn_mondec(nn2, nn_symb2),
            make_linear_pathway_nn(nn3, nn_symb3)
        )
    )
end

# Fetches the parameter ranges for the simple linear pathway model.
get_linear_pathway_ranges() = [
    :v0 => (-2.0, 1.0),
    :v => (-2.0, 1.0),
    :K => (-2.0, 1.0),
    :n => (0.0, 10.0),
    :d => (-2.0, 1.0)
]

# The model where the rate functions is fully known (the only unknown parameter is d).
function make_linear_pathway_known(params)
    # Creates the model.
    params = make_dict(params)
    v0 = params[:v0]
    v = params[:v]
    K = params[:K]
    n = params[:n]
    rs = @reaction_network known begin
        10^d, (X1, X2, X3, X4) --> 0
        10^($v0) + hillr(X1, 10^($v), 10^($K), $n), 0 --> X2
        10^($v0) + hillr(X2, 10^($v), 10^($K), $n), 0 --> X3
        10^($v0) + hillr(X3, 10^($v), 10^($K), $n), 0 --> X4
    end

    # Sets the metadata.
    MTK.get_metadata(rs)[:rates_of_interest] = [reactions(rs)[end].rate]
    MTK.get_metadata(rs)[:dependent_sps] = [[:X3]]

    return rs
end

# The model where the rate functions are parametrised. Only other unknown parameter is d.
function make_linear_pathway_params()
    # Creates the model.
    rs = @reaction_network params begin
        10^d, (X1, X2, X3, X4) --> 0
        10^v0 + hillr(X1, 10^v, 10^K, n), 0 --> X2
        10^v0 + hillr(X2, 10^v, 10^K, n), 0 --> X3
        10^v0 + hillr(X3, 10^v, 10^K, n), 0 --> X4
    end

    # Sets the metadata.
    MTK.get_metadata(rs)[:rates_of_interest] = [reactions(rs)[end].rate]
    MTK.get_metadata(rs)[:dependent_sps] = [[:X3]]

    return rs
end

# The model where the rates are fitted as a neural network. Only other unknown parameter is d.
# We have constrained the neural network to be monotonously decreasing and bounded.
function make_linear_pathway_nn_mondec_bounded(params, nn, nn_symb)
    # Creates the model.
    params = make_dict(params)
    v0 = params[:v0]
    rs = @reaction_network nn_mondec_bounded begin
        10^d, (X1, X2, X3, X4) --> 0
        10^($v0) + $nn(X1), 0 --> X2
        10^($v0) + $nn(X2), 0 --> X3
        10^($v0) + $nn(X3), 0 --> X4
    end

    # Sets the metadata.
    MTK.get_metadata(rs)[:rates_of_interest] = [reactions(rs)[end].rate]
    MTK.get_metadata(rs)[:nns] = [nn_symb]
    MTK.get_metadata(rs)[:nn_transforms] = [x -> 10^(v0) + x]
    MTK.get_metadata(rs)[:dependent_sps] = [[:X3]]

    return rs
end

# The model where the rates are fitted as a neural network. Only other unknown parameter is d.
# We have constrained the neural network to be monotonously decreasing.
function make_linear_pathway_nn_mondec(nn, nn_symb)
    # Creates the model.
    rs = @reaction_network nn_mondec begin
        10^d, (X1, X2, X3, X4) --> 0
        $nn(X1), 0 --> X2
        $nn(X2), 0 --> X3
        $nn(X3), 0 --> X4
    end

    # Sets the metadata.
    MTK.get_metadata(rs)[:rates_of_interest] = [reactions(rs)[end].rate]
    MTK.get_metadata(rs)[:nns] = [nn_symb]
    MTK.get_metadata(rs)[:dependent_sps] = [[:X3]]

    return rs
end

# The model where the rate is fitted as a neural network. Only other unknown parameter is d.
function make_linear_pathway_nn(nn, nn_symb)
    # Creates the model.
    rs = @reaction_network nn begin
        10^d, (X1, X2, X3, X4) --> 0
        $nn(X1), 0 --> X2
        $nn(X2), 0 --> X3
        $nn(X3), 0 --> X4
    end

    # Sets the metadata.
    MTK.get_metadata(rs)[:rates_of_interest] = [reactions(rs)[end].rate]
    MTK.get_metadata(rs)[:nns] = [nn_symb]
    MTK.get_metadata(rs)[:dependent_sps] = [[:X3]]

    return rs
end
