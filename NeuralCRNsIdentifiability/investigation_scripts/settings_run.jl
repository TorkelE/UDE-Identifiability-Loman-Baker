# Sets the number of processes to use for parallelisation.
num_procs = 21

# Sets the true parameter values.
ps_true = [
    :v => log10(1.1),
    :K => log10(2.0),
    :n => 3.0,
    :d => log10(0.5)
]

# Designates the model (by setting its retrieval function).
retrieval_function = get_extended_sal_constraints

# Sets the data generation settings.
u0s = [[:X => 2.0, :Y => 0.1]]
tend = 45.0
sample_n = 20
measured_sps = [:X]
σ_std = 0.3
noise_func = scaled_normal_error

# Sets the optimisation settings.
n = 200
maxiters = 20000

# Sets a seed.
seed = 12345

# Sets the profiling settings
profiling_methods = [
    CustomPLMethod(m = 50, n_starts = 40, adam_maxiters = 1000)
]
