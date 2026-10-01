# Functional and parametric identifiability for universal differential equations

Code and figures for Loman & Baker. Universal differential equations are fitted to synthetic
data generated from chemical reaction network models, and the identifiability of both the
mechanistic parameters (via profile likelihoods) and the learnt functions (via ensembles of
independent fits) is assessed.

All commands below are run from the repository root.

## Repository structure

- **`NeuralCRNsIdentifiability/`** — the Julia package. `src/` holds the model library, data
  generation, fitting, profile likelihood analysis and plotting; `investigation_scripts/`
  holds the entry points used to prepare and run scans and single runs.
- **`SingleRuns/`** — one folder per single run. A run fits every variant of one model (the
  fully specified, parameterised and UDE versions) to one synthetic dataset, and computes a
  profile likelihood for each mechanistic parameter. Each folder holds the `settings_run.jl`
  that defines it; results are written here when the run is performed.
- **`Scans/`** — one folder per scan. A scan repeats that same analysis across a range of
  settings — data quality, noise level, measured species, random seed — with one `scan_N`
  subfolder per combination. Each folder holds a `settings_run.jl` giving the shared
  settings and a `settings_scan.jl` listing the per-point overrides.
- **`Figures/`** — one folder per model, each containing the script that builds that model's
  figures, the plots it produces (`FigureComponents/`), and the assembled final versions.
- **`Project.toml`, `Manifest.toml`** — the Julia environment, pinned to the exact package
  versions used for the paper.

## Setup

Requires **Julia 1.11.5** — the version `Manifest.toml` was resolved with. Under any other
version Pkg re-resolves the manifest, replacing the pinned package versions. Use
[juliaup](https://julialang.org/install/) to manage Julia versions, and make 1.11.5 the
default so the commands below work as written:

```bash
juliaup add 1.11.5
juliaup default 1.11.5
julia --version    # must print 1.11.5
```

Install the packages once:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

`Manifest.toml` records the exact package versions used for the paper, so this reproduces
the original environment rather than resolving fresh ones. The scripts call
`Pkg.activate(".")` themselves, so `--project=.` is only needed here.

### The modified Catalyst version

The code stores information on each model in its `ReactionSystem` metadata, which requires a
mutable dictionary. That is not available in any registered Catalyst release, so the paper
used a modified version. `Manifest.toml` pins it to an archived snapshot:

```
repo-url = "https://github.com/TorkelE/Catalyst-neuralcrns-paper.git"
repo-rev = "3dad6c5cf3ce79a7801ba81e672e238e4acb3515"   # tag: neuralcrns-paper-v1
version  = "15.0.5"
```

`Pkg.instantiate()` fetches it automatically; no manual step is needed. That repository is a
frozen snapshot for this paper only — the maintained package is at
[SciML/Catalyst.jl](https://github.com/SciML/Catalyst.jl).

## Workflow

Each scan and single run ships its `settings_run.jl` (and `settings_scan.jl` for scans),
which fully specify it. Regenerating a result has three steps:

1. **Prepare** — `setup_investigation.jl` reads the stored settings and generates the
   synthetic data. For a scan it also creates one `scan_N` subfolder per scan point.
2. **Run** — `compute_run.sh` performs the fitting and profile likelihood analysis for a
   scan; `run_investigation.jl` does the same for a single run.
3. **Plot** — the scripts in `Figures/` load the completed runs and write the plots that
   were assembled into the paper figures.

The `regenerate_data` argument is what makes step 1 use the stored settings instead of the
defaults written at the top of the script.

Scans are computationally expensive.

## Figure regeneration

### Prepare the scans

```bash
julia NeuralCRNsIdentifiability/investigation_scripts/setup_investigation.jl regenerate_data extended_sal_basic_scan
julia NeuralCRNsIdentifiability/investigation_scripts/setup_investigation.jl regenerate_data extended_sal_basic_scan_alt_parameters
julia NeuralCRNsIdentifiability/investigation_scripts/setup_investigation.jl regenerate_data extended_sal_basic_scan_alt_parameters_low_conc
julia NeuralCRNsIdentifiability/investigation_scripts/setup_investigation.jl regenerate_data extended_sal_data_variation_scan
julia NeuralCRNsIdentifiability/investigation_scripts/setup_investigation.jl regenerate_data gwo_basic_scan
julia NeuralCRNsIdentifiability/investigation_scripts/setup_investigation.jl regenerate_data gwo_nn_constraints_scan
julia NeuralCRNsIdentifiability/investigation_scripts/setup_investigation.jl regenerate_data ic_feedforward_num_nn_variation_basic_scan
julia NeuralCRNsIdentifiability/investigation_scripts/setup_investigation.jl regenerate_data linear_pathway_nn_constraints_scan
julia NeuralCRNsIdentifiability/investigation_scripts/setup_investigation.jl regenerate_data mosquitoes_data_samples_scan_1
julia NeuralCRNsIdentifiability/investigation_scripts/setup_investigation.jl regenerate_data mosquitoes_data_samples_scan_2
julia NeuralCRNsIdentifiability/investigation_scripts/setup_investigation.jl regenerate_data mosquitoes_data_samples_scan_3
julia NeuralCRNsIdentifiability/investigation_scripts/setup_investigation.jl regenerate_data mosquitoes_data_samples_scan_4
julia NeuralCRNsIdentifiability/investigation_scripts/setup_investigation.jl regenerate_data mosquitoes_data_variation_scan
```

### Run the scans

Set `TAG` at the top of `investigation_scripts/compute_run.sh` to the scan you want, then:

```bash
./NeuralCRNsIdentifiability/investigation_scripts/compute_run.sh
```

To check a finished scan:

```bash
julia NeuralCRNsIdentifiability/investigation_scripts/check_scan.jl extended_sal_basic_scan
```

### Single runs

```bash
julia NeuralCRNsIdentifiability/investigation_scripts/setup_investigation.jl run regenerate_data simple_sal_basic_run
julia NeuralCRNsIdentifiability/investigation_scripts/run_investigation.jl simple_sal_basic_run 1

julia NeuralCRNsIdentifiability/investigation_scripts/setup_investigation.jl run regenerate_data extended_sal_extensive_data_run
julia NeuralCRNsIdentifiability/investigation_scripts/run_investigation.jl extended_sal_extensive_data_run 1

julia NeuralCRNsIdentifiability/investigation_scripts/setup_investigation.jl run regenerate_data modified_sir_misspecification_basic_run
julia NeuralCRNsIdentifiability/investigation_scripts/run_investigation.jl modified_sir_misspecification_basic_run 1

julia NeuralCRNsIdentifiability/investigation_scripts/setup_investigation.jl run regenerate_data modified_sir_misspecification_repeat_run
julia NeuralCRNsIdentifiability/investigation_scripts/run_investigation.jl modified_sir_misspecification_repeat_run 1

julia NeuralCRNsIdentifiability/investigation_scripts/setup_investigation.jl run regenerate_data ic_feedforward_misspecification_basic_run
julia NeuralCRNsIdentifiability/investigation_scripts/run_investigation.jl ic_feedforward_misspecification_basic_run 1
```

### Extract scan summaries

Needed by `Figures/Mosquitoes/generate_mosquitoes_figure.jl`:

```bash
julia NeuralCRNsIdentifiability/investigation_scripts/plot_scan_evaluation.jl extended_sal_data_variation_scan
julia NeuralCRNsIdentifiability/investigation_scripts/plot_scan_evaluation.jl mosquitoes_data_variation_scan
```

### Generate the figures

```bash
julia "Figures/Simple_SA_Loop/generate_simple_sa_loop_figure.jl"
julia "Figures/Extended_SA_Loop/generate_extended_sa_loop_figure.jl"
julia "Figures/Mosquitoes/generate_mosquitoes_figure.jl"
julia "Figures/IC_Feedforward/generate_ic_feedforward_figure.jl"
julia "Figures/IC_Feedforward/generate_exp_num_sup_figure.jl"
julia "Figures/Linear_Pathway/generate_linear_pathway_figure.jl"
julia "Figures/Modified_SIR/generate_modified_sir_figure.jl"
julia "Figures/Other/Optimisation_Schemes/generate_optimisation_schemes.jl"
julia "Figures/Other/Profile Likelihood Intro/generate_profile_likelihood_intro_plot.jl"
```

These write `_base` plots into each `FigureComponents/` folder; the final PDFs were
assembled from those in Inkscape (the `.svg` files alongside them).

### Which script makes which figure

Figure folders are named after the model, not after the paper figure, so several figures are
produced by a script in a folder named for a different model.

| Paper figure | File (within `Figures/`) | Generated by |
|---|---|---|
| Figure 1 | `Other/Profile Likelihood Intro/profile_likelihood_intro.pdf` | `Other/Profile Likelihood Intro/generate_profile_likelihood_intro_plot.jl` |
| Figure 2 | `Simple_SA_Loop/simple_sal_struc_nonident.pdf` | `Simple_SA_Loop/generate_simple_sa_loop_figure.jl` |
| Figure 3 | `Extended_SA_Loop/extended_sal_pract_ident.pdf` | `Extended_SA_Loop/generate_extended_sa_loop_figure.jl` |
| Figure 4 | `Mosquitoes/mosquitoes_datavar.pdf` | `Mosquitoes/generate_mosquitoes_figure.jl` |
| Figure 5 | `IC_Feedforward/ic_feedforward_2nn.pdf` | `IC_Feedforward/generate_ic_feedforward_figure.jl` |
| Figure 6 | `Linear_Pathway/ln_pathway_nn_constraints_repeat1.pdf` | `Linear_Pathway/generate_linear_pathway_figure.jl` |
| Figure 7 | `Modified_SIR/modified_sir_misspecification.pdf` | `Modified_SIR/generate_modified_sir_figure.jl` |
| Figure S1 | `Other/Optimisation_Schemes/optimiser_choice.pdf` | `Other/Optimisation_Schemes/generate_optimisation_schemes.jl` |
| Figure S2 | `Extended_SA_Loop/extended_sal_pred_uncert.pdf` | `Extended_SA_Loop/generate_extended_sa_loop_figure.jl` |
| Figure S3 | `Extended_SA_Loop/extended_sal_pract_ident_maxdata.pdf` | `Extended_SA_Loop/generate_extended_sa_loop_figure.jl` |
| Figure S4 | `Extended_SA_Loop/extended_sal_pract_ident_rep1.pdf` | `Extended_SA_Loop/generate_extended_sa_loop_figure.jl` |
| Figure S5 | `Extended_SA_Loop/extended_sal_pract_ident_rep2.pdf` | `Extended_SA_Loop/generate_extended_sa_loop_figure.jl` |
| Figure S6 | `Extended_SA_Loop/extended_sal_pract_ident_alt_ps.pdf` | `Extended_SA_Loop/generate_extended_sa_loop_figure.jl` |
| Figure S7 | `Extended_SA_Loop/gwo_pract_ident.pdf` | `Extended_SA_Loop/generate_extended_sa_loop_figure.jl` |
| Figure S8 | `Extended_SA_Loop/extended_sal_pract_ident_lowconc_ps.pdf` | `Extended_SA_Loop/generate_extended_sa_loop_figure.jl` |
| Figure S9 | `Mosquitoes/extended_sal_datavar.pdf` | `Mosquitoes/generate_mosquitoes_figure.jl` |
| Figure S10 | `Mosquitoes/extended_sal_sample_fits.pdf` | `Mosquitoes/generate_mosquitoes_figure.jl` |
| Figure S11 | `Mosquitoes/mosquitoes_sample_fits.pdf` | `Mosquitoes/generate_mosquitoes_figure.jl` |
| Figure S12 | `IC_Feedforward/ic_feedforward_exp_num_analysis.pdf` | `IC_Feedforward/generate_exp_num_sup_figure.jl` |
| Figure S13 | `IC_Feedforward/ic_feedforward_2nn_repeat.pdf` | `IC_Feedforward/generate_ic_feedforward_figure.jl` |
| Figure S14 | `Linear_Pathway/ln_pathway_nn_constraints_repeat2.pdf` | `Linear_Pathway/generate_linear_pathway_figure.jl` |
| Figure S15 | `Linear_Pathway/ln_pathway_nn_constraints_repeat3.pdf` | `Linear_Pathway/generate_linear_pathway_figure.jl` |
| Figure S16 | `Linear_Pathway/gwo_nn_constraints_repeat2.pdf` | `Linear_Pathway/generate_linear_pathway_figure.jl` |
| Figure S17 | `Linear_Pathway/gwo_nn_constraints_repeat1.pdf` | `Linear_Pathway/generate_linear_pathway_figure.jl` |
| Figure S18 | `Modified_SIR/modified_sir_misspecification_repeat.pdf` | `Modified_SIR/generate_modified_sir_figure.jl` |
| Figure S19 | `Modified_SIR/ic_feedforward_misspecification.pdf` | `Modified_SIR/generate_modified_sir_figure.jl` |

### HPC

`investigation_scripts/hpc_run.sh` is an alternative to `compute_run.sh` for SLURM clusters,
submitting a scan as an array job with one task per scan point. It can be used if more
compute is required than a single machine provides. The run tag and the array size are set at
the top of the file, and it is submitted with `sbatch`.

It loads Julia via `module load`, naming the module from the cluster used for the paper.
Point that at your own cluster's Julia 1.11.5.
