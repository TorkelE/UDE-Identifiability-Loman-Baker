# Sets the number of processes to use for parallelisation.
num_procs = 21

# Sets the true parameter values.
ps_true = [
    :v => log10(1.0),
    :K => log10(0.4),
    :n => 4.0,
    :dₘ => log10(1.0),
    :pₑ => log10(0.25),
    :dₑ => log10(0.2),
    :pₚ => log10(0.1),
    :dₚ => log10(0.05)
]

# Designates the model (by setting its retrieval function).
retrieval_function = get_gwo

# Sets the data generation settings.
u0s = [[:M => 1.0, :E => 0.1, :P => 0.1]]
tend = 30.0
sample_n = 40
measured_sps = [:P, :E]
σ_std = 0.05
noise_func = scaled_normal_error

# Sets the optimisation settings.
n = 200
maxiters = 25000

# Sets a seed.
seed = 12345

# Sets the profiling settings
profiling_methods = [
    CustomPLMethod(m = 100, n_starts = 80, adam_maxiters = 50000)
]
