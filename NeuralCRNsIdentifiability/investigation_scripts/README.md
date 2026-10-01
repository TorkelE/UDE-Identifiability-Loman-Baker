# Investigation scripts instructions

Some of this code that uses hpc (makes things faster) assumes the Oxford HPC architecture. Should be similar on many other that uses the same SLURM system.

## HPC scans: Instructions

#### Defining run settings

Modify these scripts to determine the settings of the run.
- "settings_run.jl" (for single-run, and scan based, settings).
- "settings_scan.jl" (for scan-specific settings).

(This is in teh corresponding scan/run folder in Scans/SIngleRuns)

#### Defining run character

Next, we determine the character of the run (i.e. whenever it is a scan, its name, etc.). Here, first modify ""setup_investigation.jl" with these choices. Next, on the HPC, run
```
srun -p interactive --pty /bin/bash
```
to start an interactive session, and then
```
module purge
module load Julia/1.11.3-linux-x86_64
```
to activate Julia.

Then, on hpc or on compute machines, run
```julia
julia NeuralCRNsIdentifiability/investigation_scripts/setup_investigation.jl
```
to prepare the run.

Do do stuff using the terminal only, the runtag can be given as input. You can designate that this is for a single run only using the "run" input, or to only regenerate the data using the "regenerate_data" input.
```julia
julia NeuralCRNsIdentifiability/investigation_scripts/setup_investigation.jl extended_sal_basic_scan
julia NeuralCRNsIdentifiability/investigation_scripts/setup_investigation.jl run extended_sal_basic_scan
julia NeuralCRNsIdentifiability/investigation_scripts/setup_investigation.jl regenerate_data extended_sal_basic_scan
./NeuralCRNsIdentifiability/investigation_scripts/compute_run.sh
julia NeuralCRNsIdentifiability/investigation_scripts/setup_investigation.jl run regenerate_data extended_sal_basic_scan
```

#### Running the scan (locally)

If you want to run a full scan on a compute machine, do
```
./NeuralCRNsIdentifiability/investigation_scripts/compute_run.sh
```

If you just want to run a single run (not a scan) and in Julia, instead use:
```julia
julia NeuralCRNsIdentifiability/investigation_scripts/setup_investigation.jl run regenerate_data simple_sal_basic_run
julia NeuralCRNsIdentifiability/investigation_scripts/run_investigation.jl extended_sal_basic_scan 1
```

You can omit the latest two
```julia
julia NeuralCRNsIdentifiability/investigation_scripts/run_investigation.jl
```
And they should be inferred correctly using the information in the "active_run.jls" file. 


#### Running the scan (hpc)

Submit this script to slurm (no modifications needed, will run things as an array job).
```
sbatch NeuralCRNsIdentifiability/investigation_scripts/hpc_run.sh
```
