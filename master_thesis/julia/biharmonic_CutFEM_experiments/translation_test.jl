module TranslationTest

    using Printf
    using LaTeXStrings
    using Latexify
    import Plots

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
    end

    default_size = (800, 600)
    # default_size = (400, 300)

    # Plots.pgfplotsx()
    # endfix=".tex"
    Plots.gr()
    endfix=".png"


    function sci_str(number)
        if number == 0
            return "\$ 0.0 \\cdot 10^{0} \$"
        else
            exp = round(log10(abs(number)))
            mantissa = number / 10^exp
            return @sprintf("\$%.1f \\cdot 10^{%d}\$", mantissa, exp)
        end
    end

    function translation_solve(solver; δs, L, n, γ, γg1, γg2)
        println("\nTranslation 0 to $(δs[end]);  iterations $(length(δs)); L = $L, n=$n, γ=$γ, γg1=$γg1, γg2=$γg2")

        el2s = Float64[]
        eh1s = Float64[]
        ehs_energy = Float64[]
        cond_numbers = Float64[]
        for δi in δs
            sol = solver.run_solver( n=n, u_ex=solver.u_ex, dirname=nothing,
                                    δ=δi,γ=γ, γg1=γg1, γg2=γg2, L=L)
            push!(el2s, sol.el2)
            push!(eh1s, sol.eh1)
            push!(ehs_energy, sol.eh_energy)
            push!(cond_numbers, sol.cond_number)
        end
        cond_numbers = [number > 1e23 ? 1e23 : number for number in cond_numbers] #ceiling cond numbers
        return cond_numbers, el2s, eh1s, ehs_energy
    end


    function convergence_solve(solver ;ns, γ, γg1, γg2)
        println("\n Convergence $ns, γ=$γ, γg1=$γg1, γg2=$γg2")

        el2s = Float64[]
        eh1s = Float64[]
        ehs_energy = Float64[]
        cond_numbers = Float64[]
        for ni in ns
            sol = solver.run_solver(n=ni, u_ex=solver.u_ex, dirname=nothing,
                                    γ=γ, γg1=γg1, γg2=γg2)
            push!(el2s, sol.el2)
            push!(eh1s, sol.eh1)
            push!(ehs_energy, sol.eh_energy)
            push!(cond_numbers, sol.cond_number)
        end
        cond_numbers = [number > 1e23 ? 1e23 : number for number in cond_numbers] #ceiling cond numbers
        return cond_numbers, el2s, eh1s, ehs_energy
    end

    function run_convergence(solver, param_list, dirname, prefix)
        ns = [2^3, 2^4, 2^5, 2^6, 2^7]

        # Compute results
        results = Vector{SimulationData}()
        for (params, color) in param_list
            γ, γg1, γg2 = params
            cond_numbers, el2s, eh1s, ehs_energy = convergence_solve(solver, ns=ns, γ=γ, γg1=γg1, γg2=γg2)
            push!(results, SimulationData(params, color, cond_numbers, el2s, eh1s, ehs_energy))
        end

        # Plots.gr()

        hs = 1 .// ns

        # Plot condition numbers
        p1 = Plots.plot(legend=:outertopright, size=default_size, legendtitle=L"(\gamma, \gamma_1, \gamma_2)", yscale=:log10, xscale=:log2, minorgrid=false)

        # Merge L2 error, H1 error, and Energy error into one plot
        p2 = Plots.plot(legend=:outertopright, size=default_size,legendtitle=L"(\gamma, \gamma_1, \gamma_2)", yscale=:log2, xscale=:log2, minorgrid=false)

        for sim_data in results
            γ, γg1, γg2 = sim_data.params
            label_text = L" %$(sci_str(γ)), %$(sci_str(γg1)), %$( sci_str(γg2) ) "
            Plots.plot!(p1, hs, sim_data.cond_numbers, label=label_text, color=sim_data.color)

            Plots.plot!(p2, hs, sim_data.el2s, label=label_text, color=sim_data.color, linestyle=:solid)
            Plots.plot!(p2, hs, sim_data.eh1s, label=nothing, color=sim_data.color, linestyle=:dash)
            Plots.plot!(p2, hs, sim_data.ehs_energy, label=nothing, color=sim_data.color, linestyle=:dot)

        end

        Plots.xlabel!(p1, L"$h$")
        Plots.ylabel!(p1, L"$\kappa(A)$")
        Plots.ylims!(p1, (1e5, 1e25))

        Plots.xlabel!(p2, L"$h$")
        Plots.ylabel!(p2, L"$\Vert e \Vert_{L^2,solid} $, $\Vert e \Vert_{H^1,dash} $, $\Vert e \Vert_{ah,*,dot}$")

        Plots.savefig(p1, dirname*"/$prefix"*"_cond_conv"*endfix)
        Plots.savefig(p2, dirname*"/$prefix"*"_errors_conv"*endfix)

    end


    function translation_test(solver, param_list, δs, dirname, prefix)

        function run_simulations(solver, param_list, δs, L=1.11, n=2^4)
            results = Vector{SimulationData}()

            for (params, color) in param_list
                γ, γg1, γg2 = params
                cond_numbers, el2s, eh1s, ehs_energy = translation_solve(solver, δs=δs, L=L, n=n, γ=γ, γg1=γg1, γg2=γg2)
                push!(results, SimulationData(params, color, cond_numbers, el2s, eh1s, ehs_energy))
            end
            return results
        end

        results = run_simulations(solver, param_list, δs)

        # Plot condition numbers
        p1 = Plots.plot(legend=:outertopright, size=default_size, legendtitle=L"(\gamma, \gamma_1, \gamma_2)", yscale=:log10, minorgrid=false)

        # Merge L2 error, H1 error, and Energy error into one plot
        p2 = Plots.plot(legend=:outertopright, size=default_size,legendtitle=L"(\gamma, \gamma_1, \gamma_2)", yscale=:log10, minorgrid=false)

        for sim_data in results
            γ, γg1, γg2 = sim_data.params
            label_text = L" %$(sci_str(γ)), %$(sci_str(γg1)), %$( sci_str(γg2) ) "
            Plots.plot!(p1, δs, sim_data.cond_numbers, label=label_text, color=sim_data.color)

            Plots.plot!(p2, δs, sim_data.el2s, label=label_text, color=sim_data.color, linestyle=:solid)
            Plots.plot!(p2, δs, sim_data.eh1s, label=nothing, color=sim_data.color, linestyle=:dash)
            Plots.plot!(p2, δs, sim_data.ehs_energy, label=nothing, color=sim_data.color, linestyle=:dot)

        end

        Plots.xlabel!(p1, L"$\delta$")
        Plots.ylabel!(p1, L"$\kappa(A)$")
        Plots.ylims!(p1, (1e5, 1e25))

        Plots.xlabel!(p2, L"$\delta$")
        Plots.ylabel!(p2, L"$\Vert e \Vert_{L^2,solid} $, $\Vert e \Vert_{H^1,dash} $, $\Vert e \Vert_{ah,*,dot}$")

        Plots.savefig(p1, dirname*"/$prefix"*"_cond_trans"*endfix)
        Plots.savefig(p2, dirname*"/$prefix"*"_errors_trans"*endfix)

    end

    function penalty_test(; dirname, u_ex, run_solver, iterations=5)
        solver = Solver(u_ex, run_solver)

        δ1 = 0
        L = 1.11
        n = 2^4
        h = L/n
        δ2 = sqrt(2)*h
        δs = LinRange(δ1, δ2, iterations)


        # No penlaty check
        param_list = [
            ((10., 5., 0.1), "blue"),
            ((10., 0., 0.), "red")
        ]

        prefix = "no_penalty"
        translation_test(solver, param_list, δs, dirname, prefix)
        run_convergence(solver, param_list, dirname, prefix)

    end

end
