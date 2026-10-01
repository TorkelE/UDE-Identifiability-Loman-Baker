
# Smooth transition of data of various qualities.
begin
    # Defines the data.
    n_repeats_per_setting = 40
    n_data_settings = 10
    sample_range = (60, 15)
    σ_range = (0.05, 0.40)

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
