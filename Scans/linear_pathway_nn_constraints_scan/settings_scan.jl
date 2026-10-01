# Repeats for different randomised datasets.
begin
    # Defines the scan iterator.
    scan_iterator = [
        Dict(
            :seed => 2
        ),
        Dict(
            :seed => 3,
            :profiling_methods => [CustomPLMethod(m = 40, n_starts = 40, adam_maxiters = 1000)]

        ),
        Dict(
            :seed => 7
        ),
    ]
end
