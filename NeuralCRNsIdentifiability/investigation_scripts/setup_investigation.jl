### Sets Investigation Run Settings ###
### This Section Can Be Modified ###

# General run settings (for single-runs and scans).
run_tag = "extended_sal_basic_run_alt_parameters_low_conc"
load_previous_settings = false
re_generate_data = true

# Scan-specific settings (only when a full scan is performed).
isa_scan = false
use_scanwide_data = false

# Whether to update the active run files (hpc and normal).
update_active_run_files = true

### Sets Up The Investigation Run ###
### This Section Should Not Be Modified ###

# Sets the correct folder and activates its environment.
begin
    cd("$(@__DIR__)/../..")
    using Pkg; Pkg.activate(".")
end

# If there are input arguments, we are likely on the hpc and have provided specific arguments.
if !isempty(ARGS)
    # Infers the run tag.
    run_tag = only(filter(s -> s ∉ ["run", "regenerate_data"], ARGS))

    # Infers whether this is a run or a scan.
    isa_scan = (!in("run", ARGS) && !in(run_tag, readdir("SingleRuns")))
    in("regenerate_data", ARGS) && isdir("Scans/$run_tag/scan_1") && (isa_scan = true)

    # Checks whether to regenerate the data (currently we always do this).
    re_generate_data = true

    # Checks whenever we want to load the previous settings or (or use the previous ones).
    load_previous_settings = in("regenerate_data", ARGS)

    # Checks whether to update use_scanwide_data.
    if in("regenerate_data", ARGS)
        if isfile("Scans/$run_tag/scan_1/measured_data.jls")
            use_scanwide_data = false
        elseif isfile("Scans/$run_tag/measured_data.jls")
            use_scanwide_data = true
        end
    end
end

# Prints initialisation message.
println("Preparing investigation: $run_tag")

# Fetch packages.
begin
    using NeuralCRNsIdentifiability
    import Plots: savefig
    using Serialization
end

# Initiates the run folder.
begin
    run_folder = "$(isa_scan ? "Scans" : "SingleRuns")/$run_tag"
    check_mkdir(run_folder)
    for file in readdir(run_folder)
        (file in ["settings_run.jl", "settings_scan.jl"]) || rm("$run_folder/$file", recursive = true)
    end
end

# Loads the run settings.
if load_previous_settings
    settings_dict = make_setting_dict("../../$run_folder/settings_run.jl")
    isa_scan && (scan_settings_dicts = make_scan_setting_dicts("../../$run_folder/settings_scan.jl"))
else
    settings_dict = make_setting_dict("settings_run.jl")
    isa_scan && (scan_settings_dicts = make_scan_setting_dicts("settings_scan.jl"))
end

# Saves the run settings to files.
begin
    serialize("$run_folder/settings_run.jls", settings_dict)
    load_previous_settings || cp("NeuralCRNsIdentifiability/investigation_scripts/settings_run.jl", "$run_folder/settings_run.jl"; force = true)
    if isa_scan
        serialize("$run_folder/settings_scan.jls", scan_settings_dicts)
        load_previous_settings || cp("NeuralCRNsIdentifiability/investigation_scripts/settings_scan.jl", "$run_folder/settings_scan.jl"; force = true)
    end
end

# In the case of a scan, creates all subfolder (and their overwrite settings).
if isa_scan
    for (idx, settings_override) in enumerate(scan_settings_dicts)
        isdir("$run_folder/scan_$idx/base_figures") && rm("$run_folder/scan_$idx/base_figures", recursive = true)
        check_mkdir("$run_folder/scan_$idx")
        check_mkdir("$run_folder/scan_$idx/base_figures")
        serialize("$run_folder/scan_$idx/settings_override.jls", settings_override)
    end
    use_scanwide_data && check_mkdir("$run_folder/base_figures")
else
    isdir("$run_folder/base_figures") && rm("$run_folder/base_figures", recursive=true)
    check_mkdir("$run_folder/base_figures")
end

# Generates the data. Either we should only generate a single data file. If not, we generate
# data for each part of the scan (and ensures that the common data file does not exist).
if !isa_scan || use_scanwide_data
    if re_generate_data || !isfile("$run_folder/measured_data.jls")
        md_1 = MeasuredData(settings_dict)
        serialize("$run_folder/measured_data.jls", md_1)
        savefig(plot_all_data(md_1), "$run_folder/base_figures/measured_data.png")
    end
else
    isfile("$run_folder/measured_data.jls") && rm("$run_folder/measured_data.jls")
    for (idx, settings_override) in enumerate(scan_settings_dicts)
        (idx % 5 == 0) && println("Generating data for scan $idx")
        specific_settings = combined_settings(settings_dict, settings_override)
        md_2 = MeasuredData(specific_settings)
        serialize("$run_folder/scan_$idx/measured_data.jls", md_2)
        savefig(plot_all_data(md_2), "$run_folder/scan_$idx/base_figures/measured_data.png")
    end
end

# Potentially updates the HPC file.
if update_active_run_files
    tot_scans = isdefined(Main, :scan_settings_dicts) ? length(scan_settings_dicts) : nothing
    isa_scan && make_hpc_file(run_tag, tot_scans)
    active_run_dict = Dict(:run_tag => run_tag, :isa_scan => isa_scan, :tot_scans => tot_scans)
    isa_scan && (active_run_dict[:scan_idxs] = 1)
    serialize("NeuralCRNsIdentifiability/investigation_scripts/active_run.jls", active_run_dict)
end

# Prints finalisation message.
println("Finalised preparing investigation: $run_tag")


