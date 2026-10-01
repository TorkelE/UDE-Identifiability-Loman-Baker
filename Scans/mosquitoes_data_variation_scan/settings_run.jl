# Sets the number of processes to use for parallelisation.
num_procs = 21

# Sets the true parameter values.
ps_true = [
    :r => log10(0.3),
    :d1 => log10(0.05),
    :k1 => log10(0.2),
    :l1 => log10(0.05),
    :l2 => log10(0.1),
    :l3 => log10(0.2),
    :k2 => log10(0.15),
    :d3 => log10(0.1)
]

# Designates the model (by setting its retrieval function).
retrieval_function = get_mosquitoes

# Sets the data generation settings.
u0s = [[:A => 1.0, :E => 0.0, :L => 0.0]]
tend = 80.0
sample_n = 60
measured_sps = [:A, :L]
σ_std = 0.05
noise_func = scaled_normal_error

# Sets the optimisation settings.
n = 120
maxiters = 20000

# Sets a seed.
seed = 12345

# Sets the profiling settings
profiling_methods = [
    CustomPLMethod(m = 120, n_starts = 40, adam_maxiters = 2000, p_range = (-2.0, 5.0))
]
