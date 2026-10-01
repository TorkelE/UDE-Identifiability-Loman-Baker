### Preparations ###

# Sets the correct directory. Fetches packages.
begin
    cd("$(@__DIR__)/../../..")
    using Pkg; Pkg.activate(".")
    using Plots
    import Plots.mm
end

# Set plotting defaults.
default(framestyle = :box)

### Figure generation ###

let
    x_vals = -1.0:0.001:1.0
    ident_prof(x) = -exp(abs(x)/0.1) + 1
    struct_nonident_prof(x) = 0
    pract_nonident_prof(x) = -exp((x/0.66)^4) + 1
    y_vals_ident = ident_prof.(x_vals)
    y_vals_struct_nonident = struct_nonident_prof.(x_vals)
    y_vals_pract_nonident = pract_nonident_prof.(x_vals)

    xlimit = (-1.0, 1.0); xticks = [-1.0, -0.5, 0.0, 0.5, 1.0]; ylimit = (-4.2, 0.2); yticks = [-4, -2, 0]; legend = false;
    plot(x_vals, y_vals_ident, lw = 7, color = :blue, la = 0.8)
    plt1 = plot!([x_vals[1], x_vals[end]], fill(-1.92, 2); color = :black, lw = 3.0, la = 0.5, linestyle = :dash,
        xlimit, xticks, ylimit, yticks, legend)
    plot(x_vals, y_vals_struct_nonident, lw = 7, color = :blue, la = 0.8)
    plt2 = plot!([x_vals[1], x_vals[end]], fill(-1.92, 2); color = :black, lw = 3.0, la = 0.5, linestyle = :dash,
        xlimit, xticks, ylimit, yticks, legend)
    plot(x_vals, y_vals_pract_nonident, lw = 7, color = :blue, la = 0.8)
    plt3 = plot!([x_vals[1], x_vals[end]], fill(-1.92, 2); color = :black, lw = 3.0, la = 0.5, linestyle = :dash,
        xlimit, xticks, ylimit, yticks, legend)

    lp_intro_plt = plot(plt1, plt2, plt3, size = (900, 300), layout = (1, 3), bottom_margin = 0mm, left_margin = -3mm)
    make_fig(lp_intro_plt, "Other/Profile Likelihood Intro/FigureComponents/profile_likelihood_intro_base")
end

plot(x_vals, y_vals_ident, xlimit = (-1.0, 1.0), xticks = [-1.0, -0.5, 0.0, 0.5, 1.0], ylimit = (-4.2, 0.2), yticks = [-4, -2, 0], legend = false, lw = 7, color = :blue)
