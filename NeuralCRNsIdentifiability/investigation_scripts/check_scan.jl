# Determines the scan to investigate, then runs the checking function.
begin
    run_tag = (length(ARGS) == 0) ? "extended_sal_data_variation_scan" : only(ARGS)
    cd("$(@__DIR__)/../..")
    using Pkg; Pkg.activate(".")
    using NeuralCRNsIdentifiability
    check_scan(run_tag);
    isinteractive() || exit()
end
