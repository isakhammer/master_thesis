using Plots
pgfplotsx()
default_size = (600, 200)

using CSV
using DataFrames
using LaTeXStrings
using Printf
using YAML

function main()
    println("Plotting physical CH")

    data_circle = CSV.read("physical_CH_circle/sol.csv", DataFrame)
    data_flower = CSV.read("physical_CH_flower/sol.csv", DataFrame) # Make sure to correct the file path

    p1 = plot(data_circle.ts, data_circle.e_L1_ts, size=default_size, xlabel=L"$t$", ylabel=L"$e_{L^1(\Omega)}$", legend=false, linecolor = :red)
    plot!(p1, data_flower.ts, data_flower.e_L1_ts, linecolor = :blue)

    p2 = plot(data_circle.ts[2:end], data_circle.Es[2:end],size=default_size, xscale=:log10, yscale=:log10, label = "Circle", xlabel=L"$t$", ylabel=L"$E(u_h)$", linecolor = :red)
    plot!(p2, data_flower.ts[2:end], data_flower.Es[2:end], xscale=:log10, yscale=:log10, label = "Flower", linecolor = :blue)

    p = plot(p1, p2, layout = (1, 2)) # This will arrange p1 and p2 horizontally.

    savefig(p,"physical_CH_plot.tex" )

end

main()
