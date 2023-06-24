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

    its_circle = LinRange(1, length(data_circle.ts), length(data_circle.ts))
    its_flower = LinRange(1, length(data_flower.ts), length(data_flower.ts))

    p1 = plot(its_circle, data_circle.e_L1_ts, size=default_size, xlabel=L"$t/\tau$", ylabel=L"$e_{u}$", legend=false, linecolor = :red)
    plot!(p1, its_flower, data_flower.e_L1_ts, linecolor = :blue)

    p2 = plot(its_circle[2:end], data_circle.Es[2:end],size=default_size, xscale=:log10, yscale=:log10, label = "Circle", xlabel=L"$t/\tau$", ylabel=L"$E(u_h)$", linecolor = :red)
    plot!(p2, its_flower[2:end], data_flower.Es[2:end], xscale=:log10, yscale=:log10, label = "Flower", linecolor = :blue)

    p = plot(p1, p2, layout = (1, 2)) # This will arrange p1 and p2 horizontally.

    savefig(p,"physical_CH_plot.tex" )

end

main()
