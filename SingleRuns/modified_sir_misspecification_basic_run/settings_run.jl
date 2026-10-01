# Sets the number of processes to use for parallelisation.
num_procs = 21

# Sets the true parameter values.
ps_true = [
    :β => log10(0.003),
    :γ => log10(0.005),
    :α₁ => -0.05,
    :α₂ => -0.5 # Note, we adapt due to how model is specified and with law of mass action.
]

# Designates the model (by setting its retrieval function).
retrieval_function = get_modified_sir

# Sets the data generation settings.
u0s = [[:S => 99.0, :I => 1.0, :R => 0.0]]
tend = 1000.0
sample_n = 50
measured_sps = [:I]
σ_std = 0.2
noise_func = scaled_normal_error

# Sets the optimisation settings.
n = 40
maxiters = 10000

# Sets a seed.
seed = 1

# Sets the profiling settings
profiling_methods = [
    CustomPLMethod(m = 500, n_starts = 120, adam_maxiters = 10000; p_range = (-4.0, 1.0))
]
