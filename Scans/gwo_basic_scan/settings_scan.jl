# Repeats for different randomised datasets and different data qualities.
begin
    # Defines the scan iterator.
    scan_iterator = [
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
            :seed => 2
        ),
        Dict(
            :sample_n => 30,
            :σ_std => 0.2,
            :measured_sps => [:P, :E],
            :seed => 2
        )
    ]
end
