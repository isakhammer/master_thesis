using Plots
pgfplotsx()
default_size = (600, 200)

using CSV
using DataFrames
using LaTeXStrings
using Printf
using YAML


struct Simulation
    param
    data
    dirname
    color
end

function translation_plot(sims::Vector{Simulation},dirname)

    # Plot for condition numbers
    p1 = plot(legend=:outertopright, size=default_size, legendtitle=L"(\gamma, \gamma_1, \gamma_2)",
              yscale=:log10, minorgrid=false, ymin=0.0)

    # Plot for L2 error, H1 error, and Energy error
    p2 = plot(legend=:outertopright, size=default_size,legendtitle=L"(\gamma, \gamma_1, \gamma_2)",
              yscale=:log10, minorgrid=false)

    for sim in sims
        γ, γg1, γg2 = sim.param
        label_text = L"(%$(γ), %$(γg1), %$( γg2)) "

        Plots.plot!(p1, sim.data.deltas, sim.data.cond_numbers, label=label_text, color=sim.color)
        Plots.plot!(p2, sim.data.deltas, sim.data.ehs_energy, label=nothing, color=sim.color, linestyle=:dot)
        Plots.plot!(p2, sim.data.deltas, sim.data.eh1s, label=nothing, color=sim.color, linestyle=:dash)
        Plots.plot!(p2, sim.data.deltas, sim.data.el2s, label=label_text, color=sim.color, linestyle=:solid)

    end

    scatter!(p1, [[0]], [[10^0]], color=:transparent, label=nothing, markersize=0.001)
    scatter!(p2, [[0]], [[10^-2]], color=:transparent, label=nothing, markersize=0.001)
    plot!(p2,[0], [1], linestyle = :dot, label = L"\Vert e \Vert_{a_h,*}", color = "black")
    plot!(p2,[0], [1], linestyle = :dash, label = L"\Vert e \Vert_{H^1}", color = "black")
    plot!(p2,[0], [1], linestyle = :solid, label = L"\Vert e \Vert_{L^2}", color = "black")
    plot!(p1,[0], [1], linestyle = :solid, label = L"\kappa_{\infty} (A)", color = "black")
    plot!(p2, legendfontsize=12)  # Adjust the value 12 to your desired font size
    plot!(p1, legendfontsize=12)  # Adjust the value 12 to your desired font size
    xlabel!(p1, L"\delta")
    xlabel!(p2, L"\delta")
    file1 = dirname*"/translation-cond.tex"
    file2 = dirname*"/translation-error.tex"
    println("Saved in $file1 and $file2")
    savefig(p1,file1)
    savefig(p2,file2)
end




println("Laplace Translation Test")
dirname = "translation-test/laplace_no_penalty"

# First sim
sim = "$dirname/sim-1"
data = CSV.read("$sim.csv", DataFrame)
params_yml = YAML.load_file("$sim.yml")
params = (params_yml["gamma"], params_yml["gamma1"], params_yml["gamma2"])
sim1 = Simulation(params, data, dirname, "blue")

# Second sim
sim = "$dirname/sim-2"
data = CSV.read("$sim.csv", DataFrame)
params_yml = YAML.load_file("$sim.yml")
params = (params_yml["gamma"], params_yml["gamma1"], params_yml["gamma2"])
sim2 = Simulation(params, data, dirname, "red")

# # Generate plots
sims = [sim1, sim2]
translation_plot(sims, dirname)


println("Hessian Translation Test")
dirname = "translation-test/hessian_no_penalty"

# First sim
sim = "$dirname/sim-1"
data = CSV.read("$sim.csv", DataFrame)
params_yml = YAML.load_file("$sim.yml")
params = (params_yml["gamma"], params_yml["gamma1"], params_yml["gamma2"])
sim1 = Simulation(params, data, dirname, "blue")

# Second sim
sim = "$dirname/sim-2"
data = CSV.read("$sim.csv", DataFrame)
params_yml = YAML.load_file("$sim.yml")
params = (params_yml["gamma"], params_yml["gamma1"], params_yml["gamma2"])
sim2 = Simulation(params, data, dirname, "red")

# # Generate plots
sims = [sim1, sim2]
translation_plot(sims, dirname)

