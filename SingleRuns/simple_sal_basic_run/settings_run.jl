# Sets the number of processes to use for parallelisation.
num_procs = 21

# Sets the true parameter values.
ps_true = [
    :v => log10(0.6),
    :K => log10(0.3),
    :n => 3.0,
    :d => log10(0.5)
]

# Designates the model (by setting its retrieval function).
retrieval_function = get_simple_sal

# Sets the data generation settings.
u0s = [[:X => 0.18]]
tend = 15.0
sample_n = 40
measured_sps = [:X]
σ_std = 0.05
noise_func = scaled_normal_error

# Sets the optimisation settings.
n = 100
maxiters = 10000

# Sets a seed.
seed = 1

# Sets the profiling settings
profiling_methods = [
    CustomPLMethod(m = 40, n_starts = 4, adam_maxiters = 200000)
]
