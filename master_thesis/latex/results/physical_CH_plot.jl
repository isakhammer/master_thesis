using Plots
pgfplotsx()
default_size = (500, 780)

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

    p1 = plot(its_circle[1:end-1], data_circle.Delta_uhs[2:end], size=default_size, xscale=:log10, xlabel=L"$m$", ylabel=L"$\Delta u^m$", legend=false, linecolor = :red)
    plot!(p1, its_flower[1:end-1], data_flower.Delta_uhs[2:end], xscale=:log10, yscale=:log10,  linecolor = :blue)
    # scatter!(p1, its_circle[1:end-1], data_circle.Delta_uhs[2:end], markersize=2,  markeralpha=0.5, markercolor=:red)
    # scatter!(p1, its_flower[1:end-1], data_flower.Delta_uhs[2:end], markersize=2,  markeralpha=0.5, markercolor=:blue)

    p2 = plot(its_circle[1:end-1], data_circle.delta_uhs[2:end], size=default_size, xscale=:log10, xlabel=L"$m$", ylabel=L"$\delta u^m$", legend=false, linecolor = :red)
    plot!(p2, its_flower[1:end-1], data_flower.delta_uhs[2:end], xscale=:log10, linecolor = :blue)
    # scatter!(p2, its_circle[1:end-1], data_circle.delta_uhs[2:end], markersize=2,  markeralpha=0.5, markercolor=:red)
    # scatter!(p2, its_flower[1:end-1], data_flower.delta_uhs[2:end], markersize=2,  markeralpha=0.5, markercolor=:blue)

    p3 = plot(its_circle, data_circle.Es,size=default_size, xscale=:log10,  xlabel=L"$m$", ylabel=L"$E^m$", legend=false, linecolor = :red)
    plot!(p3, its_flower, data_flower.Es, xscale=:log10,   linecolor = :blue)
    # scatter!(p3, its_circle, data_circle.Es, markersize=2,  markeralpha=0.5, markercolor=:red)
    # scatter!(p3, its_flower, data_flower.Es, markersize=2, markeralpha=0.5, markercolor=:blue)

    p4 = plot(its_circle[1:end-1], data_circle.Es[1:end-1] .- data_circle.Es[2:end], size=default_size, xscale=:log10,   xlabel=L"$m$", ylabel=L"$\delta E^m$", legend=false, linecolor = :red)
    plot!(p4, its_flower[1:end-1], data_flower.Es[1:end-1] .- data_flower.Es[2:end], xscale=:log10, yscale=:log10,       linecolor = :blue)
    # scatter!(p4, its_circle[1:end-1], data_circle.Es[1:end-1] .- data_circle.Es[2:end], markersize=2, markeralpha=0.5,markercolor=:red)
    # scatter!(p4, its_flower[1:end-1], data_flower.Es[1:end-1] .- data_flower.Es[2:end],markersize=2, markeralpha=0.5, markercolor=:blue)

    p5 = plot(its_circle,  data_circle.E1s, size=default_size,   xscale=:log10, legend=false, ylabel=L"$E_1^m$", xlabel=L"$m$", linecolor = :red)
    plot!(p5, its_flower,  data_flower.E1s, xscale=:log10,        linecolor = :blue)
    # scatter!(p5, its_circle, data_circle.E1s, markersize=2, markeralpha=0.5, markercolor=:red)
    # scatter!(p5, its_flower, data_flower.E1s, markersize=2, markeralpha=0.5, markercolor=:blue)

    p6 = plot(its_circle, data_circle.E2s, size=default_size, label = "Circle", xscale=:log10, ylabel=L"$E_2^m$", xlabel=L"$m$", linecolor = :red)
    plot!(p6, its_flower, data_flower.E2s, xscale=:log10, label = "Flower", linecolor = :blue)
    # scatter!(p6, its_circle, data_circle.E2s, markersize=2, markeralpha=0.5, markercolor=:red, label=false)
    # scatter!(p6, its_flower, data_flower.E2s, markersize=2, markeralpha=0.5, markercolor=:blue, label=false)


    p = plot(p1, p2, p3, p4, p5, p6, layout = (6, 1))
    savefig(p,"physical_CH_plot.pdf" )

end

main()
