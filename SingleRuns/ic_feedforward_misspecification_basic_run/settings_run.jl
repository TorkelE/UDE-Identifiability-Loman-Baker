# Sets the number of processes to use for parallelisation.
num_procs = 21

# Sets the true parameter values.
ps_true = [
    :d => log10(0.5),
    :v1 => log10(5.0),
    :K1 => log10(0.8),
    :n => 5,
    :v2 => log10(1.0),
    :K2 => log10(0.6)
]

# Designates the model (by setting its retrieval function).
retrieval_function = get_ic_feedforward_misspecified

# Sets the data generation settings.
u0s = [[:X => 0.0, :Y => 0.0, :Z => 0.0]]
tend = 10.0
sample_n = 25
measured_sps = [:Z, :Y]
σ_std = 0.15
noise_func = scaled_normal_error

# Sets the optimisation settings.
n = 40
maxiters = 10000

# Sets a seed.
seed = 1

# Sets the profiling settings
profiling_methods = [
    CustomPLMethod(m = 60, n_starts = 120, adam_maxiters = 10000)
]
