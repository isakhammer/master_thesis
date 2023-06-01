using Plots
pgfplotsx()
default_size = (500, 300)

using CSV
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
function compute_eoc(hs::Vector, errs1::Vector, errs2::Vector, errs3::Vector)
    eoc1 = compute_eoc(hs, errs1)
    eoc2 = compute_eoc(hs, errs2)
    eoc3 = compute_eoc(hs, errs3)
    return eoc1, eoc2, eoc3
end



function generate_plots(data,path)
    hs = 1 .// data.ns
    # Initial plot with the first data series
    p = Plots.plot(hs, data.el2s, label=L"\Vert e \Vert_{L^2}", size=default_size, legend=:outertopright, xscale=:log2, yscale=:log2, minorgrid=true)
    Plots.scatter!(p, hs, data.el2s, primary=false)

    # Add the second data series
    Plots.plot!(p, hs, data.eh1s, label=L"\Vert e \Vert_{H^1}")
    Plots.scatter!(p, hs, data.eh1s, primary=false)

    # Add the third data series
    Plots.plot!(p, hs, data.ehs_energy, label=L"\Vert e \Vert_{a_{h,*}}")
    Plots.scatter!(p, hs, data.ehs_energy, primary=false)

    # Configs
    Plots.xlabel!(p, "h")
    Plots.ylabel!(p, L"\Vert e \Vert_{}")
    Plots.plot!(p, xscale=:log2, yscale=:log2, minorgrid=true)
    Plots.plot!(p, legendfontsize=14)  # Adjust the value 12 to your desired font size

    savefig("$path-plot.tex")

    hs_str = ["1/$(n)" for n in data.ns]
    eoc_l2, eoc_eh1, eoc_eh_energy = compute_eoc(hs, data.el2s, data.eh1s, data.ehs_energy)
    header = [L"$h/L$", L"$n$", L"$\Vert e \Vert_{L^2}$", "EOC", L"$ \Vert e \Vert_{H^1}$", "EOC", L"$\Vert e \Vert_{ a_h,* }$", "EOC", L"\kappa(A)", "ndofs"]
    data = hcat(hs_str, data.ns, data.el2s,  eoc_l2, data.eh1s, eoc_eh1, data.ehs_energy, eoc_eh_energy, data.cond_numbers, data.ndofs)
    formatters = (ft_printf("%s", [1]), ft_printf("%.0f", [2]), ft_printf("%.2f", [4, 6, 8]), ft_printf("%.1E", [3, 5, 7, 9, 10]), ft_nonothing)

    filename_tex = "$path-table.tex"
    open(filename_tex, "w") do io
        pretty_table(io, data; header=header, tf = tf_latex_modern, formatters=formatters)
    end
    run(`sed -i '1d;$d' $filename_tex`) #removes \begin{table} and \end{table} env
end


# Hessian Circle
dirname = "eoc-test"
path = "$dirname/eoc-hessian-circle-L-3.11-gamma0-20-gamma1-10-gamma2-0.1/conv"
param = (20,10,0.1)
data = CSV.read("$path.csv", DataFrame)
generate_plots(data, path)


# Laplace Circle
path = "$dirname/eoc-laplace-circle-L-3.11-gamma0-20-gamma1-10-gamma2-0.1/conv"
param = (20,10,0.1)
data = CSV.read("$path.csv", DataFrame)
generate_plots(data, path)

# Laplace Flower
path = "$dirname/eoc-laplace-flower-L-3.11-gamma0-20-gamma1-10-gamma2-0.1/conv"
param = (20,10,0.1)
data = CSV.read("$path.csv", DataFrame)
generate_plots(data, path)
