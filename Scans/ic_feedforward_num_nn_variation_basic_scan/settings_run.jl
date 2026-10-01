# Sets the number of processes to use for parallelisation.
num_procs = 21

# Sets the true parameter values.
ps_true = [
    :d => log10(0.5),
    :v1 => log10(5.0),
    :K1 => log10(0.8),
    :n => 3,
    :v2 => log10(1.0),
    :K2 => log10(0.4)
]

# Designates the model (by setting its retrieval function).
retrieval_function = get_ic_feedforward

# Sets the data generation settings.
u0s = [[:X => 0.0, :Y => 0.0, :Z => 1.0]]
tend = 8.0
sample_n = 20
measured_sps = [:X, :Y, :Z]
σ_std = 0.25
noise_func = scaled_normal_error

# Sets the optimisation settings.
n = 200
maxiters = 50000

# Sets a seed.
seed = 12345

# Sets the profiling settings
profiling_methods = [
    CustomPLMethod(m = 100, n_starts = 40, adam_maxiters = 10000)
]
