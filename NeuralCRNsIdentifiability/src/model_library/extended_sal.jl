### Basic Models ###

# The extended self-activation loop. This is a model with a two components (X and Y),
# where Y activates X through some function and X activates Y linearly.

# Fetches the parameter ranges and the models for the xy self activation loop model.
function get_extended_sal(ps_true, rng;
        chain = Lux.Chain(
            Lux.Dense(1 => 3, Lux.softplus, use_bias = false),
            Lux.Dense(3 => 3, Lux.softplus, use_bias = false),
            Lux.Dense(3 => 1, Lux.softplus, use_bias = false)
        )
    )
    nn_symb, nn = make_nn(chain, rng)
    return (
        get_extended_sal_ranges(),
        (
            make_extended_sal_known(ps_true),
            make_extended_sal_params(),
            make_extended_sal_nn(nn, nn_symb)
        )
    )
end

# Fetches the parameter ranges for the xy self activation loop model.
get_extended_sal_ranges() = [
    :v => (-1.0, 1.0),
    :K => (-1.0, 1.0),
    :n => (1.0, 5.0),
    :d => (-2.0, 1.0)
]

# The model where the rate function is fully known.
function make_extended_sal_known(params)
    # Creates the model.
    params = make_dict(params)
    v = params[:v]
    K = params[:K]
    n = params[:n]
    rs = @reaction_network known begin
        hill(Y, 10^($v), 10^($K), $n), ∅ --> X
        X, 0 --> Y
        (10^d), (X,Y) --> 0
    end

    # Sets the metadata.
    MTK.get_metadata(rs)[:rates_of_interest] = [reactions(rs)[1].rate]
    MTK.get_metadata(rs)[:dependent_sps] = [[:Y]]

    return rs
end

# The model where the rate function is parametrised (with v, K, and n).
function make_extended_sal_params()
    # Creates the model.
    rs = @reaction_network params begin
        hill(Y,(10^v),(10^K),n), 0 --> X
        X, 0 --> Y
        (10^d), (X,Y) --> 0
    end

    # Sets the metadata.
    MTK.get_metadata(rs)[:rates_of_interest] = [reactions(rs)[1].rate]
    MTK.get_metadata(rs)[:dependent_sps] = [[:Y]]

    return rs
end

# The model where the rate is fitted as a neural network.
function make_extended_sal_nn(nn, nn_symb)
    # Creates the model.
    rs = @reaction_network nn begin
        $nn(Y), 0 --> X
        X, 0 --> Y
        (10^d), (X,Y) --> 0
    end

    # Sets the metadata.
    MTK.get_metadata(rs)[:rates_of_interest] = [reactions(rs)[1].rate]
    MTK.get_metadata(rs)[:nns] = [nn_symb]
    MTK.get_metadata(rs)[:dependent_sps] = [[:Y]]

    return rs
end

### Neural Network Constraints Models ###


# Fetches the parameter ranges and the models for the extended self-activation loop. Includes
# models where we have put constraints on the neural network.
function get_extended_sal_constraints(ps_true, rng;
        chain1 = Lux.Chain(
            Lux.Dense(1 => 3, Lux.softplus, use_bias = true),
            Lux.Dense(3 => 3, Lux.sigmoid_fast, use_bias = true),
            Lux.Dense(3 => 3, Lux.softplus, use_bias = true),
            Lux.Dense(3 => 1, Lux.sigmoid_fast, use_bias = true)
        ),
        chain2 = Lux.Chain(
            Lux.Dense(1 => 3, Lux.softplus, use_bias = true),
            Lux.Dense(3 => 3, Lux.sigmoid_fast, use_bias = true),
            Lux.Dense(3 => 3, Lux.softplus, use_bias = true),
            Lux.Dense(3 => 1, Lux.softplus, use_bias = true)
        ),
        chain3 = Lux.Chain(
            Lux.Dense(1 => 3, Lux.softplus, use_bias = false),
            Lux.Dense(3 => 3, Lux.softplus, use_bias = false),
            Lux.Dense(3 => 1, Lux.softplus, use_bias = false)
        ),
    )
    nn_symb1, nn1 = make_nn(chain1, rng)
    nn_symb2, nn2 = make_nn(chain2, rng)
    nn_symb3, nn3 = make_nn(chain3, rng)
    return (
        get_extended_sal_ranges(),
        (
            make_extended_sal_known(ps_true),
            make_extended_sal_params(),
            make_extended_sal_nn_moninc_bounded(ps_true, nn1, nn_symb1),
            make_extended_sal_nn_moninc(nn2, nn_symb2),
            make_extended_sal_nn(nn3, nn_symb3)
        )
    )
end

# The model where the rates are fitted as a neural network. Only other unknown parameter is d.
# We have constrained the neural network to be monotonously decreasing and bounded.
function make_extended_sal_nn_moninc_bounded(params, nn, nn_symb)
    # Creates the model.
    params = make_dict(params)
    v = params[:v]

    # Creates the model.
    rs = @reaction_network nn_moninc_bounded begin
        (10^($v)) * ($nn(Y)), 0 --> X
        X, 0 --> Y
        (10^d), (X,Y) --> 0
    end

    # Sets the metadata.
    MTK.get_metadata(rs)[:rates_of_interest] = [reactions(rs)[1].rate]
    MTK.get_metadata(rs)[:nns] = [nn_symb]
    MTK.get_metadata(rs)[:nn_transforms] = [y -> (10^(v)) * y]
    MTK.get_metadata(rs)[:dependent_sps] = [[:Y]]

    return rs
end

# The model where the rates are fitted as a neural network. Only other unknown parameter is d.
# We have constrained the neural network to be monotonously decreasing.
function make_extended_sal_nn_moninc(nn, nn_symb)
    # Creates the model.
    rs = @reaction_network nn_moninc begin
        $nn(Y), 0 --> X
        X, 0 --> Y
        (10^d), (X,Y) --> 0
    end

    # Sets the metadata.
    MTK.get_metadata(rs)[:rates_of_interest] = [reactions(rs)[1].rate]
    MTK.get_metadata(rs)[:nns] = [nn_symb]
    MTK.get_metadata(rs)[:dependent_sps] = [[:Y]]

    return rs
end
