# Sets the number of processes to use for parallelisation.
num_procs = 41

# Sets the true parameter values.
ps_true = [
    :v0 => log10(0.05),
    :v => log10(1.00),
    :K => log10(0.4),
    :n => 7.0,
    :d => log10(1.0)
]

# Designates the model (by setting its retrieval function).
retrieval_function = get_linear_pathway

# Sets the data generation settings.
u0s = [[:X1 => 1.0, :X2 => 1.0, :X3 => 1.0, :X4 => 1.0]]
tend = 5.0
sample_n = 50
measured_sps = [:X4]
σ_std = 0.15
noise_func = scaled_normal_error

# Sets the optimisation settings.
n = 200
maxiters = 50000

# Sets a seed.
seed = 12345

# Sets the profiling settings
profiling_methods = [
    CustomPLMethod(m = 30, n_starts = 30, adam_maxiters = 1000)
]
