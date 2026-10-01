# A model with a single component (X) which activates its own production according to some function.

# Fetches the parameter ranges and the models for the simple self activation loop model.
function get_simple_sal(ps_true, rng;
        chain = Lux.Chain(
            Lux.Dense(1 => 3, Lux.softplus, use_bias = false),
            Lux.Dense(3 => 3, Lux.softplus, use_bias = false),
            Lux.Dense(3 => 1, Lux.softplus, use_bias = false)
        )
    )
    nn_symb, nn = make_nn(chain, rng)
    return (
        get_simple_sal_ranges(),
        (
            make_simple_sal_known(ps_true),
            make_simple_sal_params(),
            make_simple_sal_nn(nn, nn_symb)
        )
    )
end

# Fetches the parameter ranges for the simple self activation loop model.
get_simple_sal_ranges() = [
    :v => (-1.0, 1.0),
    :K => (-1.0, 1.0),
    :n => (1.0, 5.0),
    :d => (-2.0, 3.0)
]

# The model where the rate function is fully known.
function make_simple_sal_known(params)
    # Creates the model.
    params = make_dict(params)
    v = params[:v]
    K = params[:K]
    n = params[:n]
    rs = @reaction_network known begin
        hill(X, 10^($v), 10^($K), $n), ∅ --> X
        10^d, X --> ∅
    end

    # Sets the metadata.
    MTK.get_metadata(rs)[:rates_of_interest] = [reactions(rs)[1].rate]
    MTK.get_metadata(rs)[:dependent_sps] = [[:X]]

    return rs
end

# The model where the rate function is parametrised (with K and n).
function make_simple_sal_params()
    # Creates the model.
    rs = @reaction_network params begin
        hill(X, 10^v, 10^K, n), ∅ --> X
        10^d, X --> ∅
    end

    # Sets the metadata.
    MTK.get_metadata(rs)[:rates_of_interest] = [reactions(rs)[1].rate]
    MTK.get_metadata(rs)[:dependent_sps] = [[:X]]

    return rs
end

# The model where the rate is fitted as a neural network.
function make_simple_sal_nn(nn, nn_symb)
    # Creates the model.
    rs = @reaction_network nn begin
        $nn(X), ∅ --> X
        10^d, X --> ∅
    end

    # Sets the metadata.
    MTK.get_metadata(rs)[:rates_of_interest] = [reactions(rs)[1].rate]
    MTK.get_metadata(rs)[:nns] = [nn_symb]
    MTK.get_metadata(rs)[:dependent_sps] = [[:X]]

    return rs
end
