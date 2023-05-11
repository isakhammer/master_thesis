module TranslationTest

    using LaTeXStrings
    using Latexify
    import Plots

    # default_size = (800, 600)
    default_size = (400, 300)

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

    function translation_solve(;δs, L, n, γ, γg1, γg2)
        println("\nTranslation $δ1 to $δ2;  iterations $(length(δs)); L = $L, n=$n, γ=$γ, γg1=$γg1, γg2=$γg2")

        el2s = Float64[]
        eh1s = Float64[]
        ehs_energy = Float64[]
        cond_numbers = Float64[]
        for δi in δs
            sol = Solver.run( n=n, u_ex=u_ex, dirname=dirname,
                             δ=δi,γ=γ, γg1=γg1, γg2=γg2, L=L)
            push!(el2s, sol.el2)
            push!(eh1s, sol.eh1)
            push!(ehs_energy, sol.eh_energy)
            push!(cond_numbers, sol.cond_number)
        end
        cond_numbers = [number > 1e23 ? 1e23 : number for number in cond_numbers] #ceiling cond numbers
        return cond_numbers, el2s, eh1s, ehs_energy
    end

    function run_simulations(param_list, δs, L, n)
        results = Vector{SimulationData}()

        for (params, color) in param_list
            γ, γg1, γg2 = params
            cond_numbers, el2s, eh1s, ehs_energy = translation_solve(δs=δs, L=L, n=n, γ=γ, γg1=γg1, γg2=γg2)
            push!(results, SimulationData(params, color, cond_numbers, el2s, eh1s, ehs_energy))
        end

        return results
    end

    function convergence_solve(;ns, γ, γg1, γg2)
        println("\n Convergence $ns, n=$n, γ=$γ, γg1=$γg1, γg2=$γg2")

        el2s = Float64[]
        eh1s = Float64[]
        ehs_energy = Float64[]
        cond_numbers = Float64[]
        for ni in ns
            sol = Solver.run( n=ni, u_ex=u_ex,
                             γ=γ, γg1=γg1, γg2=γg2)
            push!(el2s, sol.el2)
            push!(eh1s, sol.eh1)
            push!(ehs_energy, sol.eh_energy)
            push!(cond_numbers, sol.cond_number)
        end
        cond_numbers = [number > 1e23 ? 1e23 : number for number in cond_numbers] #ceiling cond numbers
        return cond_numbers, el2s, eh1s, ehs_energy
    end

    function run_convergence(param_list, dirname, prefix)
        ns = [2^3, 2^4, 2^5, 2^6, 2^7]

        # Compute results
        results = Vector{SimulationData}()
        for (params, color) in param_list
            γ, γg1, γg2 = params
            cond_numbers, el2s, eh1s, ehs_energy = convergence_solve(ns=ns, γ=γ, γg1=γg1, γg2=γg2)
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


    function create_plot_from_results(results, δs, dirname, prefix)


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

        Plots.pgfplotsx()
        Plots.savefig(p1, dirname*"/$prefix"*"_cond_trans"*endfix)
        Plots.savefig(p2, dirname*"/$prefix"*"_errors_trans"*endfix)

    end

    function translation_test(; dirname, u_ex )
        iterations = 1000
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

        results = run_simulations(param_list, δs, L, n)

        prefix = "no_penalty"
        create_plot_from_results(results, δs, dirname, prefix)
        run_convergence(param_list, dirname, prefix)


    end

    function parameter_sweep()
        # Base params
        σ = [10^6, 10^5, 10^4, 10^3, 10^2, 10^-2, 10^-3, 10^-4, 10^-5, 10^-6]

        # Parameter Sweep γg2
        param_list = [
                      ((10., 5., 0.1 * σ[1]), "blue"),
                      ((10., 5., 0.1 * σ[2]), "red"),
                      ((10., 5., 0.1 * σ[3]), "green"),
                      ((10., 5., 0.1 * σ[4]), "purple"),
                      ((10., 5., 0.1 * σ[5]), "orange"),
                      ((10., 5., 0.1 * σ[6]), "cyan"),
                      ((10., 5., 0.1 * σ[7]), "brown"),
                      ((10., 5., 0.1 * σ[8]), "magenta"),
                      ((10., 5., 0.1 * σ[9]), "pink"),
                      ((10., 5., 0.1 * σ[10]), "lime"),
                     ]

        results = run_simulations(param_list, δs, L, n)
        prefix = "sweep_g2"
        create_plot_from_results(results, δs, dirname, prefix)
        run_convergence(param_list, dirname, prefix)

        # Parameter Sweep γg1
        param_list = [
                      ((10., 5. * σ[1] , 0.1), "blue"),
                      ((10., 5. * σ[2] , 0.1), "red"),
                      ((10., 5. * σ[3] , 0.1), "green"),
                      ((10., 5. * σ[4] , 0.1), "purple"),
                      ((10., 5. * σ[5] , 0.1), "orange"),
                      ((10., 5. * σ[6] , 0.1), "cyan"),
                      ((10., 5. * σ[7] , 0.1), "brown"),
                      ((10., 5. * σ[8] , 0.1), "magenta"),
                      ((10., 5. * σ[9] , 0.1), "pink"),
                      ((10., 5. * σ[10], 0.1), "lime"),
                     ]

        results = run_simulations(param_list, δs, L, n)
        prefix = "sweep_g1"
        create_plot_from_results(results, δs, dirname, prefix)
        run_convergence(param_list, dirname, prefix)

        # Parameter Sweep γg1 and γg2
        param_list = [
                      ((10., 5.0*σ[1] , 0.1 * σ[1] ), "blue"),
                      ((10., 5.0*σ[2] , 0.1 * σ[2] ), "red"),
                      ((10., 5.0*σ[3] , 0.1 * σ[3] ), "green"),
                      ((10., 5.0*σ[4] , 0.1 * σ[4] ), "purple"),
                      ((10., 5.0*σ[5] , 0.1 * σ[5] ), "orange"),
                      ((10., 5.0*σ[6] , 0.1 * σ[6] ), "cyan"),
                      ((10., 5.0*σ[7] , 0.1 * σ[7] ), "brown"),
                      ((10., 5.0*σ[8] , 0.1 * σ[8] ), "magenta"),
                      ((10., 5.0*σ[9] , 0.1 * σ[9] ), "pink"),
                      ((10., 5.0*σ[10], 0.1 * σ[10]), "lime"),
                     ]

        results = run_simulations(param_list, δs, L, n)
        prefix = "sweep_g1_g2"
        create_plot_from_results(results, δs, dirname, prefix)
        run_convergence(param_list, dirname, prefix)

    end

end
