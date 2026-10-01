# A incoherent feedforward loop model. The production parameter for X has been set to 1 to
# achieve structural identifiability. The two hill functions uses the same n.

# Fetches the parameter ranges and the models for the incoherent feedforward loop model.
function get_ic_feedforward(ps_true, rng;
        chain = Lux.Chain(
            Lux.Dense(1 => 3, Lux.sigmoid_fast, use_bias = false),
            Lux.Dense(3 => 3, Lux.softplus, use_bias = false),
            Lux.Dense(3 => 1, Lux.softplus, use_bias = false)
        )
    )
    nn_symb1, nn1 = make_nn(chain, rng; nn_idx = 1)
    nn_symb2, nn2 = make_nn(chain, rng; nn_idx = 2)
    return (
        get_ic_feedforward_ranges(),
        (
            make_ic_feedforward_known(ps_true),
            make_ic_feedforward_params_params(),
            make_ic_feedforward_nn_params(nn1, nn_symb1),
            make_ic_feedforward_params_nn(nn2, nn_symb2),
            make_ic_feedforward_nn_nn(nn1, nn2, nn_symb1, nn_symb2)
        )
    )
end

# Fetches the parameter ranges for the simple self activation loop model.
get_ic_feedforward_ranges() = [
    :d => (-2.0, 1.0),
    :v1 => (-2.0, 1.0),
    :K1 => (-2.0, 1.0),
    :n => (1.0, 5.0),
    :v2 => (-2.0, 1.0),
    :K2 => (-2.0, 1.0)
]

# The model where the rate functions are fully known (the only unknown parameter is d).
function make_ic_feedforward_known(params)
    # Creates the model.
    params = make_dict(params)
    d = params[:d]
    v1 = params[:v1]
    K1 = params[:K1]
    n = params[:n]
    v2 = params[:v2]
    K2 = params[:K2]
    rs = @reaction_network known begin
        1.0, 0 --> X
        hillr(X, 10^($v1), 10^($K1), $n), 0 --> Y
        hill(X, 10^($v2), 10^($K2), $n), Y --> Z
        10^d, (X,Y,Z) --> ∅
    end

    # Sets the metadata.
    MTK.get_metadata(rs)[:rates_of_interest] = [
        reactions(rs)[2].rate,
        reactions(rs)[3].rate,
    ]
    MTK.get_metadata(rs)[:dependent_sps] = [[:X], [:X]]
    MTK.get_metadata(rs)[:is_nn_rate] = [false, false]

    return rs
end

# The model where both rate functions are parametrised. Only other unknown parameter is d.
function make_ic_feedforward_params_params()
    # Creates the model.
    rs = @reaction_network params_params begin
        1.0, 0 --> X
        hillr(X, 10^v1, 10^K1, n), 0 --> Y
        hill(X, 10^v2, 10^K2, n), Y --> Z
        10^d, (X,Y,Z) --> ∅
    end

    # Sets the metadata.
    MTK.get_metadata(rs)[:rates_of_interest] = [
        reactions(rs)[2].rate,
        reactions(rs)[3].rate,
    ]
    MTK.get_metadata(rs)[:dependent_sps] = [[:X], [:X]]
    MTK.get_metadata(rs)[:is_nn_rate] = [false, false]

    return rs
end

# The model where the first rate functions is parametrised and the second is fitted as a neural network.
# Only other unknown parameter is d.
function make_ic_feedforward_params_nn(nn2, nn_symb2)
    # Creates the model.
    rs = @reaction_network params_nn begin
        1.0, 0 --> X
        hillr(X, 10^v1, 10^K1, n), 0 --> Y
        $nn2(X), Y --> Z
        10^d, (X,Y,Z) --> ∅
    end

    # Sets the metadata.
    MTK.get_metadata(rs)[:rates_of_interest] = [
        reactions(rs)[2].rate,
        reactions(rs)[3].rate,
    ]
    MTK.get_metadata(rs)[:nns] = [nn_symb2]
    MTK.get_metadata(rs)[:dependent_sps] = [[:X], [:X]]
    MTK.get_metadata(rs)[:is_nn_rate] = [false, true]

    return rs
end

# The model where the first rate functions is fitted as a neural network and the second is parametrised.
# Only other unknown parameter is d.
function make_ic_feedforward_nn_params(nn1, nn_symb1)
    # Creates the model.
    rs = @reaction_network nn_params begin
        1.0, 0 --> X
        $nn1(X), 0 --> Y
        hill(X, 10^v2, 10^K2, n), Y --> Z
        10^d, (X,Y,Z) --> ∅
    end

    # Sets the metadata.
    MTK.get_metadata(rs)[:rates_of_interest] = [
        reactions(rs)[2].rate,
        reactions(rs)[3].rate,
    ]
    MTK.get_metadata(rs)[:nns] = [nn_symb1]
    MTK.get_metadata(rs)[:dependent_sps] = [[:X], [:X]]
    MTK.get_metadata(rs)[:is_nn_rate] = [true, false]

    return rs
end

# The model where both rate functions are fitted as neural networks. Only other unknown parameter is d.
function make_ic_feedforward_nn_nn(nn1, nn2, nn_symb1, nn_symb2)
    # Creates the model.
    rs = @reaction_network nn_nn begin
        1.0, 0 --> X
        $nn1(X), 0 --> Y
        $nn2(X), Y --> Z
        10^d, (X,Y,Z) --> ∅
    end

    # Sets the metadata.
    MTK.get_metadata(rs)[:rates_of_interest] = [
        reactions(rs)[2].rate,
        reactions(rs)[3].rate,
    ]
    MTK.get_metadata(rs)[:nns] = [nn_symb1, nn_symb2]
    MTK.get_metadata(rs)[:dependent_sps] = [[:X], [:X]]
    MTK.get_metadata(rs)[:is_nn_rate] = [true, true]
    return rs
end

### the Case With Misspecified Model ###

# Fetches the parameter ranges and the models for the modified SIR model.
function get_ic_feedforward_misspecified(ps_true, rng;
        chain = Lux.Chain(
            Lux.Dense(1 => 3, Lux.softplus, use_bias = false),
            Lux.Dense(3 => 3, Lux.softplus, use_bias = false),
            Lux.Dense(3 => 1, Lux.softplus, use_bias = false)
        )
    )
    nn_symb, nn = make_nn(chain, rng)
    return (
        get_ic_feedforward_ranges(),
        (
            make_ic_feedforward_known(ps_true),
            make_ic_feedforward_params_params(),
            make_ic_feedforward_nn_params(nn, nn_symb),
            make_ic_feedforward_wrong()
        )
    )
end

# The model where both rate functions are parametrised. Only other unknown parameter is d.# The model where the rate function is a simple parameter β (yielding the wrong model). The
# only unknown parameter is γ.
function make_ic_feedforward_wrong()
    # Creates the model.
    rs = @reaction_network wrong begin
        1.0, 0 --> X
        mmr(X, 10^v1, 10^K1), 0 --> Y
        hill(X, 10^v2, 10^K2, n), Y --> Z
        10^d, (X,Y,Z) --> ∅
    end

    # Sets the metadata.
    MTK.get_metadata(rs)[:rates_of_interest] = [
        reactions(rs)[2].rate,
        reactions(rs)[3].rate,
    ]
    MTK.get_metadata(rs)[:dependent_sps] = [[:X], [:X]]
    MTK.get_metadata(rs)[:is_nn_rate] = [false, false]

    return rs
end
