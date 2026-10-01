# Smooth transition of data of various qualities.
begin
    # Defines the data.
    n_repeats_per_setting = 40
    n_data_settings = 10
    sample_range = (100, 30)
    σ_range = (0.01, 0.07)

    # Loops through all cases and adds them to the scan iterator.
    scan_iterator = []
    sample_vals = round.(Int64, collect(range(sample_range[1], sample_range[2], length = n_data_settings)))
    σ_vals = range(σ_range[1], σ_range[2], length = n_data_settings)
    for (sample_n, σ_std) in zip(sample_vals, σ_vals)
        for i in 1:n_repeats_per_setting
            push!(scan_iterator, Dict(:sample_n => sample_n, :σ_std => σ_std, :seed => i))
        end
    end
end

# We cannot do all.
scan_iterator = [scan_iterator[1:6]; scan_iterator[121:126]; scan_iterator[241:246]; scan_iterator[361:366]]
