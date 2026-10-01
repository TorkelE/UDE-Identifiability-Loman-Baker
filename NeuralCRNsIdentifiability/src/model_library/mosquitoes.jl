# A mosquito model. Two achieve full structural identifiability one parameter must be set to 1,
# however, in practise, to achieve practical identifiability additional parameters are fixed.

# Fetches the parameter ranges and the models for the mosquito model.
function get_mosquitoes(ps_true, rng;
        chain = Lux.Chain(
            Lux.Dense(1 => 2, Lux.softplus),
            Lux.Dense(2 => 2, Lux.softplus),
            Lux.Dense(2 => 1, Lux.softplus)
        )
    )
    nn_symb, nn = make_nn(chain, rng)
    return (
        get_mosquitoes_ranges(),
        (
            make_mosquitoes_known(ps_true),
            make_mosquitoes_params(ps_true),
            make_mosquitoes_nn(ps_true, nn, nn_symb)
        )
    )
end

# Fetches the parameter ranges for the mosquito model.
get_mosquitoes_ranges() = [
    :r => (-2.0, 1.0),
    :d1 => (-2.0, 1.0),
    :d3 => (-2.0, 1.0),
    :k1 => (-2.0, 1.0),
    :k2 => (-2.0, 1.0),
    :l1 => (-2.0, 1.0),
    :l2 => (-2.0, 1.0),
    :l3 => (-2.0, 1.0)
]

# The model where the rate function is fully known (the only unknown parameter is r).
function make_mosquitoes_known(params)
    # Creates the model.
    params = make_dict(params)
    d1 = params[:d1]
    k1 = params[:k1]
    l1 = params[:l1]
    l2 = params[:l2]
    l3 = params[:l3]
    k2 = params[:k2]
    d3 = params[:d3]
    rs = @reaction_network known begin
        10^(r), A --> A + E
        10^($d1), E --> ∅
        10^($k1), E --> L
        (10^($l1)) + (10^($l2))*L + (10^($l3))*L^2, L --> 0
        10^($k2), L --> A
        10^($d3), A --> 0
    end

    # Sets the metadata.
    MTK.get_metadata(rs)[:rates_of_interest] = [reactions(rs)[4].rate]
    MTK.get_metadata(rs)[:dependent_sps] = [[:L]]

    return rs
end

# The model where the rate function is parametrised (with l1, l2, and l3). Only other unknown parameter is r.
function make_mosquitoes_params(params)
    # Creates the model.
    params = make_dict(params)
    d1 = params[:d1]
    k1 = params[:k1]
    k2 = params[:k2]
    d3 = params[:d3]
    rs = @reaction_network params begin
        10^(r), A --> A + E
        10^($d1), E --> ∅
        10^($k1), E --> L
        (10^(l1)) + (10^(l2))*L + (10^(l3))*L^2, L --> 0
        10^($k2), L --> A
        10^($d3), A --> 0
    end

    # Sets the metadata.
    MTK.get_metadata(rs)[:rates_of_interest] = [reactions(rs)[4].rate]
    MTK.get_metadata(rs)[:dependent_sps] = [[:L]]

    return rs
end

# The model where the rate is fitted as a neural network. Only other unknown parameter is r.
function make_mosquitoes_nn(params, nn, nn_symb)
    # Creates the model.
    params = make_dict(params)
    d1 = params[:d1]
    k1 = params[:k1]
    k2 = params[:k2]
    d3 = params[:d3]
    rs = @reaction_network nn begin
        10^(r), A --> A + E
        10^($d1), E --> ∅
        10^($k1), E --> L
        $nn(L), L --> 0
        10^($k2), L --> A
        10^($d3), A --> 0
    end

    # Sets the metadata.
    MTK.get_metadata(rs)[:rates_of_interest] = [reactions(rs)[4].rate]
    MTK.get_metadata(rs)[:nns] = [nn_symb]
    MTK.get_metadata(rs)[:dependent_sps] = [[:L]]

    return rs
end
