using Plots
pgfplotsx()
default_size = (600, 600)

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

    println( data_circle.Delta_uhs)

    p1 = plot(its_circle[1:end-1], data_circle.Delta_uhs[2:end], size=default_size, yscale=:log10, xscale=:log10, xlabel=L"$t/\tau$", ylabel=L"$\Delta u^m$", legend=false, linecolor = :red)
    plot!(p1, its_flower[1:end-1], data_flower.Delta_uhs[2:end], xscale=:log10,  linecolor = :blue)

    p2 = plot(its_circle[1:end-1], data_circle.delta_uhs[2:end], size=default_size, xscale=:log10, xlabel=L"$t/\tau$", ylabel=L"$\delta u^m$", legend=false, linecolor = :red)
    plot!(p2, its_flower[1:end-1], data_flower.delta_uhs[2:end], xscale=:log10, linecolor = :blue)

    p3 = plot(its_circle, data_circle.Es,size=default_size, xscale=:log10, yscale=:log10,  xlabel=L"$t/\tau$", ylabel=L"$E^m$", legend=false, linecolor = :red)
    plot!(p3, its_flower, data_flower.Es, xscale=:log10, yscale=:log10,  linecolor = :blue)

    p4 = plot(its_circle[1:end-1], data_circle.Es[1:end-1] .- data_circle.Es[2:end], size=default_size, xscale=:log10,  label = "Circle", xlabel=L"$t/\tau$", ylabel=L"$\delta E^m$", linecolor = :red)
    plot!(p4, its_flower[1:end-1], data_flower.Es[1:end-1] .- data_flower.Es[2:end], xscale=:log10, yscale=:log10, label = "Flower", linecolor = :blue)

    p = plot(p1, p2, p3, p4, layout = (4, 1))

    savefig(p,"physical_CH_plot.tex" )

end

main()
