using Plots
pgfplotsx()
default_size = (500, 300)

using CSV
using DataFrames
using LaTeXStrings
using Printf




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

    savefig("$path-eoc.tex")
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
