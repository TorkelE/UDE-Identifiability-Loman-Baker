# Repeats for different randomised datasets and different data qualities.
begin
    scan_iterator = [
        # Dataset 1 (this is run with more compute in likelihood profiling, as this dataset required a more fine grained grid to look good).
        Dict(
            :sample_n => 20,
            :σ_std => 0.1,
            :measured_sps => [:X],
            :seed => 1,
            :profiling_methods => [CustomPLMethod(m = 50, n_starts = 150, adam_maxiters = 20000)]
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

        # Dataset 2.
        Dict(
            :sample_n => 20,
            :σ_std => 0.1,
            :measured_sps => [:X],
            :seed => 222,
            :profiling_methods => [CustomPLMethod(m = 50, n_starts = 150, adam_maxiters = 20000)]
        ),
        Dict(
            :sample_n => 20,
            :σ_std => 0.1,
            :measured_sps => [:X, :Y],
            :seed => 222,
            :profiling_methods => [CustomPLMethod(m = 50, n_starts = 150, adam_maxiters = 20000)]
        ),
        Dict(
            :sample_n => 20,
            :σ_std => 0.3,
            :measured_sps => [:X, :Y],
            :seed => 222,
            :profiling_methods => [CustomPLMethod(m = 40, n_starts = 40, adam_maxiters = 5000)]
        ),

        # Dataset 3.
        Dict(
            :sample_n => 20,
            :σ_std => 0.1,
            :measured_sps => [:X],
            :seed => 333,
            :profiling_methods => [CustomPLMethod(m = 40, n_starts = 40, adam_maxiters = 5000)]
        ),
        Dict(
            :sample_n => 20,
            :σ_std => 0.1,
            :measured_sps => [:X, :Y],
            :seed => 333,
            :profiling_methods => [CustomPLMethod(m = 40, n_starts = 40, adam_maxiters = 5000)]
        ),
        Dict(
            :sample_n => 20,
            :σ_std => 0.3,
            :measured_sps => [:X, :Y],
            :seed => 333,
            :profiling_methods => [CustomPLMethod(m = 40, n_starts = 40, adam_maxiters = 5000)]
        ),
    ]
end
