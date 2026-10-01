module NeuralCRNsIdentifiability

### Fetch Required Packages ###

# Model creation, simulation, and neural networks.
using Catalyst
using Lux
using OrdinaryDiffEq
import ComponentArrays
import ComponentArrays: ComponentArray
import ModelingToolkit as MTK
import ModelingToolkitNeuralNets: SymbolicNeuralNetwork

# Parameter fitting and profile likelihood analysis
using Optimization
using OptimizationNLopt
using QuasiMonteCarlo
using LikelihoodProfiler
using LineSearches
import Optimization: LBFGS
import OptimizationNLopt
import OptimizationOptimJL
import OptimizationOptimisers: Adam

# Miscellaneous.
using CSV
using Distributions
using KernelDensity
using Plots
using Random
using Serialization
using StableRNGs

### Run Files & Export Functions ###

# A library of implemented models.
export get_simple_sal, get_extended_sal, get_extended_sal_constraints, get_mosquitoes, get_ic_feedforward, get_ic_feedforward_misspecified, get_linear_pathway, get_gwo, get_gwo_nn, get_gwo_constraints, get_modified_sir
include("model_library/simple_sal.jl")
include("model_library/extended_sal.jl")
include("model_library/mosquitoes.jl")
include("model_library/ic_feedforward.jl")
include("model_library/linear_pathway.jl")
include("model_library/gwo.jl")
include("model_library/modified_sir.jl")

# Functions for generating synthetic data.
export MeasuredData, import_data
export normal_error
include("data_generation.jl")

# Functions for fitting parameters to data.
export make_eoptprob, make_callback, fit_model
export summarise_optimisation
include("parameter_fitting.jl")

# Functions for performing profile likelihood-based identifiability analysis.
export compute_lprofiles, compute_lprofile
export OptimisationPLMethod, IntegrationPLMethod, CustomPLMethod
include("profile_likelihood_analysis.jl")

# For investigating the results of runs and scans on compute machines and HPCs.
# Has to be run before plotting, as plotting uses some structures declared here.
export load_investigation, make_summarisation_plot, plot_combined_lprofiles, plot_all_combined_lprofiles
export plot_combined_fitted_funcs, plot_fitted_sims_combined, plot_condensed_info
export plot_combined_nn_arch_eval, plot_nn_arch_fit_eval, plot_nn_arch_overfit_eval
export plot_full_data_quality_eval, plot_profile_likelihood_intervals, plot_nonmeasured_mean_err, plot_mean_func_diff
export plot_condensed_info_progression
include("saved_run_investigation.jl")

# Functions for plotting various results.
export plot_data, plot_all_data
export plot_fitted_simulation, plot_all_fitted_simulation, plot_loss_traces
export plot_fitted_func, plot_fitted_funcs, plot_fitted_func!, plot_fitted_funcs!
export plot_lprofile, plot_lprofile!, plot_lprofiles
export plotting_colors_ic, plotting_colors_linpathway
export extract_datavar_scan_evaluation, extract_nn_architecture_grid_data
include("plotting.jl")

# Various utility functions
export is_nn_model, num_nns, get_dependent_sps, getname, is_valid, make_fig
include("utility.jl")

# Functionality related to scanning stuff on HPC (or via scripts).
export check_scan
export run_investigation, make_setting_dict, make_scan_setting_dicts, interpret_input, make_hpc_file
export check_mkdir, combined_settings, retrieve_models
include("hpc_scan_support.jl")

end # module NeuralCRNsIdentifiability
