# Repeats for different randomised datasets and different data qualities.
begin
    scan_iterator = [
        # Seed 1.
        Dict(
            :sample_n => 20,
            :σ_std => 0.1,
            :measured_sps => [:X],
            :seed => 1
        ),
        Dict(
            :sample_n => 20,
            :σ_std => 0.1,
            :measured_sps => [:X, :Y],
            :seed => 1
        ),
        Dict(
            :sample_n => 20,
            :σ_std => 0.3,
            :measured_sps => [:X, :Y],
            :seed => 1
        ),
    ]
end
