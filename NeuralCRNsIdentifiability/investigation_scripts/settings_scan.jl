# Select and enable a single one of these.

# Different data qualities.
if true
    # Defines the scan iterator.
    scan_iterator = [
        # Seed 1.
        Dict(
            :sample_n => 30,
            :σ_std => 0.05,
            :measured_sps => [:P],
            :seed => 1
        ),
        Dict(
            :sample_n => 30,
            :σ_std => 0.05,
            :measured_sps => [:P, :E],
            :seed => 1
        ),
        Dict(
            :sample_n => 30,
            :σ_std => 0.2,
            :measured_sps => [:P, :E],
            :seed => 1
        ),

        # Seed 2.
        Dict(
            :sample_n => 30,
            :σ_std => 0.05,
            :measured_sps => [:P],
            :seed => 2
        ),
        Dict(
            :sample_n => 30,
            :σ_std => 0.05,
            :measured_sps => [:P, :E],
            :seed => 2
        ),
        Dict(
            :sample_n => 30,
            :σ_std => 0.2,
            :measured_sps => [:P, :E],
            :seed => 2
        ),

        # Seed 3.
        Dict(
            :sample_n => 30,
            :σ_std => 0.05,
            :measured_sps => [:P],
            :seed => 3
        ),
        Dict(
            :sample_n => 30,
            :σ_std => 0.05,
            :measured_sps => [:P, :E],
            :seed => 3
        ),
        Dict(
            :sample_n => 30,
            :σ_std => 0.2,
            :measured_sps => [:P, :E],
            :seed => 3
        ),
    ]
end


# Grids of neural network architectures.
if false
    # Overwriting error check.
    (@isdefined scan_iterator) && error("The scan iterator is about to be overwritten.")

    # Designates the layers to try.
    n_repeats_per_setting = 5
    n_inputs = 1
    n_outputs = 1
    n_layers = 1:8
    layer_width = 1:6

    # Loops through all alternatives and adds them to the scan iterator.
    scan_iterator = []
    for n in n_layers, w in layer_width
        for _ in 1:n_repeats_per_setting
            internal_layers = [
                Lux.Dense(w => w, Lux.softplus, use_bias = false)
                for i in 1:n
            ]
            push!(scan_iterator, Dict(:chain => Lux.Chain(
                    Lux.Dense(n_inputs => w, Lux.softplus, use_bias = false),
                    internal_layers...,
                    Lux.Dense(w => n_outputs, Lux.sigmoid_fast, use_bias = false)
                )
            ))
        end
    end
end

# Smooth transition of data of various qualities.
if false
    # Overwriting error check.
    (@isdefined scan_iterator) && error("The scan iterator is about to be overwritten.")

    # Defines the data.
    n_repeats_per_setting = 10
    n_data_settings = 10
    sample_range = (80, 20)
    σ_range = (0.01, 0.40)

    # Loops through all cases and adds them to the scan iterator.
    scan_iterator = []
    sample_vals = round.(Int64, collect(range(sample_range[1], sample_range[2], length = n_data_settings)))
    σ_vals = logrange(σ_range[1], σ_range[2], length = n_data_settings)
    for (sample_n, σ_std) in zip(sample_vals, σ_vals)
        for _ in 1:n_repeats_per_setting
            push!(scan_iterator, Dict(:sample_n => sample_n, :σ_std => σ_std, :seed => rand(1:10000)))
        end
    end
end

# Overwriting error check.
(@isdefined scan_iterator) || error("The scan iterator has not been defined.")
