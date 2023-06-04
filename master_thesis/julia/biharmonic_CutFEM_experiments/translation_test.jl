include("biharmonic_CutCIP_laplace.jl")
include("biharmonic_CutCIP_hessian.jl")

using Dates
using Test

module TranslationTest

    using Printf
    import Plots
    using DataFrames
    using YAML
    using CSV
    using LaTeXStrings
    using Latexify
    import Gridap

    Plots.gr()
    endfix = ".png"

    default_size = (800, 600)
    # default_size = (400, 300)
    struct Solver
        u_ex
        run_solver
    end

    struct SimulationData
        params::Tuple{Float64, Float64, Float64}
        color::String
        cond_numbers::Vector{Float64}
        el2s::Vector{Float64}
        eh1s::Vector{Float64}
        ehs_energy::Vector{Float64}
        graphics::Vector{Any}
    end

    function sci_str(number)
        # Converts float to string of scientific notation
        if number == 0
            return "\$ 0.0 \\cdot 10^{0} \$"
        else
            exp = round(log10(abs(number)))
            mantissa = number / 10^exp
            return @sprintf("\$%.1f \\cdot 10^{%d}\$", mantissa, exp)
        end
    end

    function translation_solve(solver; δs, L, n, γ, γg1, γg2)
        println("Translation 0 to $(δs[end]);  iterations $(length(δs)); L = $L, n=$n, γ=$γ, γg1=$γg1, γg2=$γg2")

        el2s = Float64[]
        eh1s = Float64[]
        ehs_energy = Float64[]
        cond_numbers = Float64[]
        graphics = Vector{Any}()
        N = length(δs)

        # Run simulation while translating the grid
        for i in 1:N
            if i%10 == 0
                println("Iteration $i/$N" )
            end
            δi = δs[i]
            sol, graphic = solver.run_solver( n=n, u_ex=solver.u_ex, dirname=nothing,
                                            δ=δi,γ=γ, γg1=γg1, γg2=γg2, L=L)
            push!(el2s, sol.el2)
            push!(eh1s, sol.eh1)
            push!(ehs_energy, sol.eh_energy)
            push!(cond_numbers, sol.cond_number)
            push!(graphics, graphic)
        end

        cond_numbers = [number > 1e23 ? 1e23 : number for number in cond_numbers] #ceiling cond numbers
        return cond_numbers, el2s, eh1s, ehs_energy, graphics
    end


    function translation_test(solver, param_list, δs, L, n, dirname)

        # run simulations
        results = Vector{SimulationData}()
        for (params, color) in param_list
            γ, γg1, γg2 = params
            cond_numbers, el2s, eh1s, ehs_energy, graphics = translation_solve(solver, δs=δs, L=L, n=n, γ=γ, γg1=γg1, γg2=γg2)
            push!(results, SimulationData(params, color, cond_numbers, el2s, eh1s, ehs_energy, graphics))
        end

        # Plot for condition numbers
        p1 = Plots.plot(legend=:outertopright, size=default_size, legendtitle=L"(\gamma, \gamma_1, \gamma_2)", yscale=:log10, minorgrid=false)

        # Merge L2 error, H1 error, and Energy error into one plot
        p2 = Plots.plot(legend=:outertopright, size=default_size,legendtitle=L"(\gamma, \gamma_1, \gamma_2)", yscale=:log10, minorgrid=false)

        for i in 1:length(results)
            sim_data = results[i]

            γ, γg1, γg2 = sim_data.params
            cond_max = maximum(sim_data.cond_numbers)
            el2_max = maximum(sim_data.el2s)
            eh1_max = maximum(sim_data.eh1s)
            eh_energy_max = maximum(sim_data.ehs_energy)

            println("γ, γg1, γg2 = $γ, $γg1, $γg2:  Max (cond, el2, eh1, eh1s, ehs_energy) = $cond_max, $el2_max, $eh1_max, $eh_energy_max")

            label_text = L" %$(sci_str(γ)), %$(sci_str(γg1)), %$( sci_str(γg2) ) "
            Plots.plot!(p1, δs, sim_data.cond_numbers, label=label_text, color=sim_data.color)
            Plots.plot!(p2, δs, sim_data.el2s, label=label_text, color=sim_data.color, linestyle=:solid)
            Plots.plot!(p2, δs, sim_data.eh1s, label=nothing, color=sim_data.color, linestyle=:dash)
            Plots.plot!(p2, δs, sim_data.ehs_energy, label=nothing, color=sim_data.color, linestyle=:dot)

            parameters = Dict("gamma" => γ,
                              "gamma1" => γg1,
                              "gamma2" => γg2,
                              "L" => L,
                              "n" => n,
                              "iterations" => length(δs)
                             )
            YAML.write_file("$dirname/sim-$i.yml", parameters)
            CSV.write(dirname*"/sim-$i.csv", DataFrame(deltas = δs, cond_numbers = sim_data.cond_numbers, el2s = sim_data.el2s, eh1s = sim_data.eh1s, ehs_energy = sim_data.ehs_energy), delim=',')

            # Plotting moving grid (boundary is standstill)
            vtkdirname = "$dirname/graphics"
            mkpath(vtkdirname)
            Gridap.writevtk(sim_data.graphics[1].Γ, vtkdirname*"/sim-$i-boundary.vtu")
            Gridap.createpvd(vtkdirname*"/sim-$i") do pvd
                N = length(δs)
                for j in 1:N
                    δi = δs[i]
                    pvd[i] = Gridap.createvtk(sim_data.graphics[i].Ω_bg, vtkdirname*"/sim-$i-delta-$j.vtu")
                end
            end

        end

        Plots.xlabel!(p1, L"$\delta$")
        Plots.ylabel!(p1, L"$\kappa(A)$")
        Plots.ylims!(p1, (1e5, 1e25))

        Plots.xlabel!(p2, L"$\delta$")
        Plots.ylabel!(p2, L"$\Vert e \Vert_{L^2,solid} $, $\Vert e \Vert_{H^1,dash} $, $\Vert e \Vert_{ah,*,dot}$")

        Plots.savefig(p1, dirname*"/cond-trans.png")
        Plots.savefig(p2, dirname*"/errors-trans.png")

        return results
    end


end

function main()
    l, m, r = (1, 1, 1)
    u_ex(x) = (x[1]^2 + x[2]^2 - 1)^2*sin(m*( 2π/l )*x[1])*cos(r*( 2π/l )*x[2])

    # Parameters
    iterations = 500
    δ1 = 0
    L = 2.7
    n = 2^4
    h = L/n
    δ2 = 2*sqrt(2)*h # two squares
    γ, γg1, γg2 = (20, 10, 0.5)
    δs = LinRange(δ1, δ2, iterations)


    # No penalty comparison
    function run_penalty_test(solver_function, dirname)
        param_list = [
                      ((γ, γg1, γg2), "blue"),
                      ((γ, 0, 0), "red")
                     ]

        # Construct solver
        solver = TranslationTest.Solver(u_ex, solver_function)
        results = TranslationTest.translation_test(solver, param_list, δs, L, n, dirname)
        sim_data_ghost_penalty, sim_data_no_penalty = results

    end

    maindir = "figures/translation-test"
    if isdir(maindir)
        rm(maindir; recursive=true)
        mkdir(maindir)
    end

    dirname = "$maindir/laplace_no_penalty"
    mkpath(dirname)
    run_penalty_test(SolverLaplace.run, dirname)

    # Make figure env
    dirname = "$maindir/hessian_no_penalty"
    mkpath(dirname)
    run_penalty_test(SolverHessian.run, dirname)

end


@time main()

