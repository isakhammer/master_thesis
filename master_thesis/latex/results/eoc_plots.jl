using Plots
pgfplotsx()
default_size = (300, 300)

using CSV
using YAML
using DataFrames
using PrettyTables
using LaTeXStrings
using Printf



# Function to compute EOC and prepend nothing for each pair of hs and errs
function compute_eoc(hs::Vector, errs::Vector)
    eoc = log.(errs[1:end-1] ./ errs[2:end]) ./ log.(hs[1:end-1] ./ hs[2:end])
    return [NaN; eoc]
end

# Function to compute EOC for three pairs of hs and errs vectors
function compute_eoc(hs::Vector, errs1::Vector, errs2::Vector, errs3::Vector, errs4::Vector)
    eoc1 = compute_eoc(hs, errs1)
    eoc2 = compute_eoc(hs, errs2)
    eoc3 = compute_eoc(hs, errs3)
    eoc4 = compute_eoc(hs, errs4)
    return eoc1, eoc2, eoc3, eoc4
end


function generate_plots(data, path)
    hs = 2 ./ data.ns

    # Set up an empty plot with the desired properties
    p = Plots.plot(size=default_size, legend=:outertopright, xscale=:log2, yscale=:log2, minorgrid=true)

    # Add the energy data series
    Plots.plot!(p, hs, data.ehs_energy, label=L"\Vert e \Vert_{a_{h,*}}")
    Plots.scatter!(p, hs, data.ehs_energy, primary=false)

    # Add the L2 data series
    Plots.plot!(p, hs, data.el2s, label=L"\Vert e \Vert_{L^2}")
    Plots.scatter!(p, hs, data.el2s, primary=false)

    # Add the H1 data series
    Plots.plot!(p, hs, data.eh1s, label=L"\Vert e \Vert_{H^1}")
    Plots.scatter!(p, hs, data.eh1s, primary=false)

    # Add the third data series
    hs_hat =  1 ./ ( 2 .^(1:11) )
    println(hs_hat)
    Plots.plot!(p, hs_hat, 300*hs_hat.^1, color="grey", label=L"O(h)")
    Plots.plot!(p, hs_hat, 100*hs_hat.^2, color="black", label=L"O(h^2)")

    # Configs
    Plots.xlabel!(p, "h")
    Plots.plot!(p, minorgrid=false)
    Plots.plot!(p, legendfontsize=12)  # Adjust the value 12 to your desired font size

    # Create a new plot for condition numbers
    p_cond = Plots.plot(size=default_size, legend=:outertopright, xscale=:log2, yscale=:log10, minorgrid=true)
    Plots.plot!(p_cond, hs, data.cond_numbers, color="red", label=L"\kappa_{\infty}(A)")
    Plots.scatter!(p_cond, hs, data.cond_numbers, color="red", primary=false)

    # Add O(h^-4) series into the condition numbers plot
    Plots.plot!(p_cond, hs_hat, 10^5*hs_hat.^(-4), color="black", label=L"O(h^{-4})")

    Plots.xlabel!(p_cond, "h")
    Plots.plot!(p_cond, minorgrid=false)
    Plots.plot!(p_cond, legendfontsize=12)  # Adjust the value 12 to your desired font size

    # Combine the two plots into a single plot with a vertical layout
    p_combined = Plots.plot(p, p_cond, layout=(2, 1), size=default_size)

    savefig(p_combined, "$path-plot.tex")

    hs_str = ["1/$(n)" for n in data.ns]
    eoc_l2, eoc_eh1, eoc_eh_energy, eoc_cond_numbers = compute_eoc(hs, data.el2s, data.eh1s, data.ehs_energy, data.cond_numbers)
    header = [L"$h/L$", L"$n$", L"$\Vert e \Vert_{L^2}$", "EOC", L"$ \Vert e \Vert_{H^1}$", "EOC", L"$\Vert e \Vert_{ a_h,* }$", "EOC", L"\kappa_{\infty}(A)", "EOC", "ndofs"]
    data = hcat(hs_str, data.ns, data.el2s, eoc_l2, data.eh1s, eoc_eh1, data.ehs_energy, eoc_eh_energy, data.cond_numbers, eoc_cond_numbers, data.ndofs)
    formatters = (ft_printf("%s", [1]), ft_printf("%.0f", [2]), ft_printf("%.2f", [4, 6, 8, 10]), ft_printf("%.1E", [3, 5, 7, 9, 11]), ft_nonothing)

    filename_tex = "$path-table.tex"
    open(filename_tex, "w") do io
        pretty_table(io, data; header=header, tf = tf_latex_modern, formatters=formatters)
    end
    run(`sed -i '1d;$d' $filename_tex`) #removes \begin{table} and \end{table} env

end


dirname = "eoc-test"
params = YAML.load_file("$dirname/parameters.yml")
println(params)

# Laplace Flower
path = "$dirname/eoc-laplace-flower/conv"
data = CSV.read("$path.csv", DataFrame)
generate_plots(data, path)

# Laplace Circle
path = "$dirname/eoc-laplace-circle/conv"
data = CSV.read("$path.csv", DataFrame)
generate_plots(data, path)

# Hessian Circle
path = "$dirname/eoc-hessian-circle/conv"
data = CSV.read("$path.csv", DataFrame)
generate_plots(data, path)

