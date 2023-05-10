using Dates
using Printf
import Plots

using LaTeXStrings
using Latexify
using PrettyTables

module Solver
    using Gridap
    using LinearAlgebra
    using PROPACK
    using GridapEmbedded
    using GridapGmsh
    using Parameters
    import Gridap: ∇

    # α(x) = x[1]^2 + x[2]^2
    α(x) = 1

    function man_sol(u_ex)
        ∇u_ex(x) = ∇(u_ex)(x)
        ∇Δu_ex(x) = ∇(Δ(u_ex))(x)
        f(x) = Δ(Δ(u_ex))(x)+ α(x)⋅u_ex(x)
        return u_ex, f, ∇u_ex, ∇Δu_ex
    end

    @with_kw struct Solution
        el2
        eh1
        eh_energy
        cond_number
        ndof
    end


    function run(; n, u_ex, dirname, L=1.11, δ=0.0, γ=10, γg1=5, γg2=0.1)

        order = 2
        u_ex, f, ∇u_ex, ∇Δu_ex = man_sol(u_ex)

        # Background model (translated)
        θ_δ = 0*( π/4 )
        r_δ = δ*Point(cos(θ_δ),sin(θ_δ))
        pmin = Point(-L, -L) + r_δ
        pmax = Point(L , L) + r_δ
        bgorigin = ( pmin + pmax )/2

        R  = 1.0
        println("Sim: order=$order, n=$n, bg (L,L)=($(round(L, digits=2)),$(round(L, digits=2))), bgorigin=($(round(bgorigin[1], digits=2)),$(round(bgorigin[2], digits=2))), disk R=$(round(R, digits=1))")
        # Background model
        partition = (n,n)
        bgmodel = CartesianDiscreteModel(pmin, pmax, partition)
        geo = disk(R)

        # Cut the background model
        cutgeo = cut(bgmodel, geo)
        cutgeo_facets = cut_facets(bgmodel,geo)

        # Set up interpolation mesh and function spaces
        Ω_act = Triangulation(cutgeo, ACTIVE)

        # Construct function spaces
        V = TestFESpace(Ω_act, ReferenceFE(lagrangian, Float64, order), conformity=:H1)
        U = TrialFESpace(V)

        # Set up integration meshes, measures and normals
        Ω = Triangulation(cutgeo, PHYSICAL)
        Γ = EmbeddedBoundary(cutgeo)
        Λ = SkeletonTriangulation(cutgeo_facets)
        Fg = GhostSkeleton(cutgeo)

        # Set up integration measures
        degree = 2*order
        dΩ   = Measure(Ω, degree)
        dΓ   = Measure(Γ, degree)
        dΛ   = Measure(Λ, degree) # F_int
        dFg  = Measure(Fg, degree)

        # Set up normal vectors
        n_Γ = get_normal_vector(Γ)
        n_Λ = get_normal_vector(Λ)
        n_Fg = get_normal_vector(Fg)

        # Mesh size
        h = L/n

        function mean_n(u,n)
            return 0.5*( u.plus⋅n.plus + u.minus⋅n.minus )
        end

        function mean_nn(u,n)
            return 0.5*( n.plus⋅ ∇∇(u).plus⋅ n.plus + n.minus ⋅ ∇∇(u).minus ⋅ n.minus )
        end

        function jump_nn(u,n)
            return ( n.plus⋅ ∇∇(u).plus⋅ n.plus - n.minus ⋅ ∇∇(u).minus ⋅ n.minus )
        end

        # Inner facets
        # a(u,v) =( ∫( ∇∇(v)⊙∇∇(u) + α⋅(v⊙u) )dΩ
        #          + ∫(-mean_nn(v,n_Λ)⊙jump(∇(u)⋅n_Λ) - mean_nn(u,n_Λ)⊙jump(∇(v)⋅n_Λ))dΛ
        #          + ∫(-( n_Γ ⋅ ∇∇(v)⋅ n_Γ )⊙∇(u)⋅n_Γ - ( n_Γ ⋅ ∇∇(u)⋅ n_Γ )⊙∇(v)⋅n_Γ)dΓ
        #          + ∫((γ/h)⋅jump(∇(u)⋅n_Λ)⊙jump(∇(v)⋅n_Λ))dΛ + ∫((γ/h)⋅ ∇(u)⊙n_Γ⋅∇(v)⊙n_Γ )dΓ
        #         )

        a_CIP(u,v) = ∫(u*v)*dΩ + ( ∫(Δ(v)⊙Δ(u))dΩ
                    + ∫(-mean(Δ(v))⊙jump(∇(u)⋅n_Λ) - mean(Δ(u))⊙jump(∇(v)⋅n_Λ) + (γ/h)⋅jump(∇(u)⋅n_Λ)⊙jump(∇(v)⋅n_Λ))dΛ
                    + ∫(-Δ(v)⊙∇(u)⋅n_Γ - Δ(u)⊙∇(v)⋅n_Γ + (γ/h)⋅ ∇(u)⊙n_Γ⋅∇(v)⊙n_Γ )dΓ
                )

        # Define linear form
        # Notation: g_1 = ∇u_ex⋅n_Γ, g_2 = ∇Δu_ex⋅n_Γ
        g_1 = ∇u_ex⋅n_Γ
        g_2 = ∇Δu_ex⋅n_Γ
        l(v) = (∫( f*v ) * dΩ
                +  ∫(-(g_2⋅v))dΓ
                + ∫(g_1⊙(-Δ(v) + (γ/h)*∇(v)⋅n_Γ)) * dΓ
               )

        g(u,v) = h^(-2)*( ∫( (γg1*h)*jump(n_Fg⋅∇(u))*jump(n_Fg⋅∇(v)) ) * dFg +
                         ∫( (γg2*h^3)*jump_nn(u,n_Fg)*jump_nn(v,n_Fg) ) * dFg)

        A(u,v) = a_CIP(u,v) + g(u,v)


        # Assemble of system
        op = AffineFEOperator(A, l, U, V)
        uh = solve(op)
        A_mat =  get_matrix(op)
        ndof = size(A_mat)[1]
        cond_number = ( 1/sqrt(ndof) )*cond(A_mat,Inf)

        u_inter = interpolate(u_ex, V) # remove?
        e = u_inter - uh
        el2 = sqrt(sum( ∫(e*e)dΩ ))

        # TODO: Add α into ∫(e⊙e)*dΩ
        eh_energy = sqrt(sum( ∫(e⊙e)dΩ + ∫( ∇∇(e)⊙∇∇(e) )dΩ
                      + ( γ/h ) * ∫(jump(∇(e)⋅n_Λ) ⊙ jump(∇(e)⋅n_Λ))dΛ
                      + ( h/γ ) * ∫(mean_nn(e,n_Λ) ⊙ mean_nn(e,n_Λ))dΛ
                      + ( γ/h ) * ∫((∇(e)⋅n_Γ) ⊙ (∇(e)⋅n_Γ))dΓ
                      + ( h/γ ) * ∫(( n_Γ ⋅ ∇∇(e)⋅ n_Γ ) ⊙ ( n_Γ ⋅ ∇∇(e)⋅ n_Γ ))dΓ
                     ))

        eh1 = sqrt(sum( ∫( e⊙e + ∇(e)⊙∇(e) )dΩ ))

        vtkdirname =dirname*"/g_$(γ)_g1_$(γg1)_g2_$(γg2)_order_$(order)_n_$n"
        mkpath(vtkdirname)

        # Write out models and computational domains for inspection
        writevtk(bgmodel,   vtkdirname*"/model")
        writevtk(Ω,         vtkdirname*"/Omega")
        writevtk(Ω_act,     vtkdirname*"/Omega_act")
        writevtk(Λ,         vtkdirname*"/Lambda")
        writevtk(Γ,         vtkdirname*"/Gamma")
        writevtk(Fg,        vtkdirname*"/Fg")
        writevtk(Λ,         vtkdirname*"/jumps",      cellfields=["jump_u"=>jump(uh)])
        writevtk(Ω,         vtkdirname*"/sol",        cellfields=["e"=>e, "uh"=>uh, "u"=>u_inter])

        sol = Solution( el2=el2, eh1=eh1, eh_energy=eh_energy,
                        cond_number=cond_number, ndof=ndof)
        return sol
    end

end # module


function generate_figures(;ns, el2s, eh1s, ehs_energy, cond_numbers, ndofs, dirname::String)
    filename = dirname*"/conv"

    hs = 1 .// ns
    # hs_str =  latexify.(hs)

    compute_eoc(hs, errs) = log.(errs[1:end-1]./errs[2:end])./log.(hs[1:end-1]./hs[2:end])
    eoc_l2 = compute_eoc(hs, el2s)
    eoc_eh1 = compute_eoc(hs,eh1s)
    eoc_eh_energy = compute_eoc(hs,ehs_energy)

    eoc_l2 =  [nothing; eoc_l2]
    eoc_eh1 =  [nothing; eoc_eh1]
    eoc_eh_energy =  [nothing; eoc_eh_energy]
    header = [L"$n$", L"$\Vert e \Vert_{L^2}$", "EOC", L"$ \Vert e \Vert_{H^1}$", "EOC", L"$\Vert e \Vert_{ a_h,* }$", "EOC", "Cond number", "ndofs"]
    data = hcat(ns, el2s,  eoc_l2, eh1s, eoc_eh1, ehs_energy, eoc_eh_energy, cond_numbers, ndofs)

    formatters = (ft_printf("%.0f", [1]), ft_printf("%.2f", [3, 5, 7]), ft_printf("%.1E", [2, 4, 6, 8, 9]), ft_nonothing)

    open(filename*".tex", "w") do io
        pretty_table(io, data, header=header, backend=Val(:latex ), formatters = formatters )
    end
    minimal_header = ["n", "L2", "EOC", "H1", "EOC", "a_h", "EOC", "cond", "const", "ndofs"]
    data = hcat(ns, el2s,  eoc_l2, eh1s, eoc_eh1, ehs_energy, eoc_eh_energy, cond_numbers, cond_numbers.*hs.^4, ndofs)
    formatters = ( ft_printf("%.0f",[1,10]), ft_printf("%.2f",[3,5,7]), ft_printf("%.1E",[2,4,6,8,9]), ft_nonothing )
    pretty_table(data, header=minimal_header, formatters =formatters )

    open(filename*".txt", "w") do io
        pretty_table(io, data, header=minimal_header, backend=:text, formatters=formatters)
    end
    # Initial plot with the first data series
    p = Plots.plot(hs, el2s, label=L"\Vert e \Vert_{L^2}", legend=:bottomright, xscale=:log2, yscale=:log2, minorgrid=true)
    Plots.scatter!(p, hs, el2s, primary=false)

    # Add the second data series
    Plots.plot!(p, hs, eh1s, label=L"\Vert e \Vert_{H^1}")
    Plots.scatter!(p, hs, eh1s, primary=false)

    # Add the third data series
    Plots.plot!(p, hs, ehs_energy, label=L"\Vert e \Vert_{a_{h,*}}")
    Plots.scatter!(p, hs, ehs_energy, primary=false)

    # Configs
    Plots.xlabel!(p, "h")
    Plots.plot!(p, xscale=:log2, yscale=:log2, minorgrid=true)
    Plots.plot!(p, legendfontsize=14)  # Adjust the value 12 to your desired font size


    # Save the plot as a .png file using the GR backend
    Plots.gr()
    # Plots.pgfplotsx()
    Plots.savefig(p, filename*"_plot.png")
    # Plots.savefig(p, filename*"_plot.tex")
end

function convergence_analysis(; ns, dirname, u_ex)

    el2s = Float64[]
    eh1s = Float64[]
    ehs_energy = Float64[]
    cond_numbers = Float64[]
    ndofs = Float64[]

    println("Convergence test", ns)
    for n in ns

        sol = Solver.run(n=n, u_ex=u_ex, dirname=dirname)

        push!(el2s, sol.el2)
        push!(eh1s, sol.eh1)
        push!(ehs_energy, sol.eh_energy)
        push!(cond_numbers, sol.cond_number)
        push!(ndofs, sol.ndof)
    end
    generate_figures(ns=ns, el2s=el2s, eh1s=eh1s, ehs_energy=ehs_energy,
                     cond_numbers=cond_numbers, ndofs=ndofs, dirname=dirname)
end


# Custom struct to hold parameters, color, and data
struct SimulationData
    params::Tuple{Float64, Float64, Float64}
    color::String
    cond_numbers::Vector{Float64}
    el2s::Vector{Float64}
    eh1s::Vector{Float64}
    ehs_energy::Vector{Float64}
end

function translation_test(; dirname, u_ex )
    iterations = 5
    δ1 = 0
    δ2 = 0.2
    L = 1.11 + δ2
    n=2^6
    δs = LinRange(δ1, δ2, iterations)


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

    function create_plot_from_results(results, δs, dirname)
        function sci_str(number)
            if number == 0
                return "\$ 0.0 \\cdot 10^{0} \$"
            else
                exp = round(log10(abs(number)))
                mantissa = number / 10^exp
                return @sprintf("\$%.1f \\cdot 10^{%d}\$", mantissa, exp)
            end
        end

        Plots.gr()

        # Plot condition numbers
        p1 = Plots.plot(legend=:outertopright, legendtitle=L"(\gamma, \gamma_1, \gamma_2)", yscale=:log10, minorgrid=false)

        # Plot L2 error
        p2 = Plots.plot(legend=:outertopright, legendtitle=L"(\gamma, \gamma_1, \gamma_2)", yscale=:log10, minorgrid=false)

        # Plot H1 error
        p3 = Plots.plot(legend=:outertopright, legendtitle=L"(\gamma, \gamma_1, \gamma_2)", yscale=:log10, minorgrid=false)

        # Plot Energy error
        p4 = Plots.plot(legend=:outertopright, legendtitle=L"(\gamma, \gamma_1, \gamma_2)", yscale=:log10, minorgrid=false)

        for sim_data in results
            γ, γg1, γg2 = sim_data.params
            label_text = L" %$(sci_str(γ)), %$(sci_str(γg1)), %$( sci_str(γg2) ) "
            Plots.plot!(p1, δs, sim_data.cond_numbers, label=label_text, color=sim_data.color)
            Plots.scatter!(p1, δs, sim_data.cond_numbers, primary=false, markerstrokealpha=0.4, markersize=3, color=sim_data.color)

            Plots.plot!(p2, δs, sim_data.el2s, label=label_text, color=sim_data.color)
            Plots.scatter!(p2, δs, sim_data.el2s, primary=false, markerstrokealpha=0.4, markersize=3, color=sim_data.color)

            Plots.plot!(p3, δs, sim_data.eh1s, label=label_text, color=sim_data.color)
            Plots.scatter!(p3, δs, sim_data.eh1s, primary=false, markerstrokealpha=0.4, markersize=3, color=sim_data.color)

            Plots.plot!(p4, δs, sim_data.ehs_energy, label=label_text, color=sim_data.color)
            Plots.scatter!(p4, δs, sim_data.ehs_energy, primary=false, markerstrokealpha=0.4, markersize=3, color=sim_data.color)
        end

        Plots.xlabel!(p1, L"$\delta$")
        Plots.ylabel!(p1, L"$\kappa(A)$")
        Plots.ylims!(p1, (1e5, 1e25))

        Plots.xlabel!(p2, L"$\delta$")
        Plots.ylabel!(p2, L"el2")

        Plots.xlabel!(p3, L"$\delta$")
        Plots.ylabel!(p3, L"eh1")

        Plots.xlabel!(p4, L"$\delta$")
        Plots.ylabel!(p4, L"e_ah")

        Plots.savefig(p1, dirname*"/cond_trans.png")
        Plots.savefig(p2, dirname*"/l2_trans.png")
        Plots.savefig(p3, dirname*"/h1_trans.png")
        Plots.savefig(p4, dirname*"/ah_trans.png")
        return
    end

    # No penlaty check
    param_list = [
        ((10., 5., 0.1), "blue"),
        ((10., 0., 0.), "red")
    ]

    plotdir = dirname*"/no_penalty"
    mkpath(plotdir)
    results = run_simulations(param_list, δs, L, n)
    create_plot_from_results(results, δs, plotdir)


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

    plotdir = dirname*"/parameter_sweep_g2"
    mkpath(plotdir)
    results = run_simulations(param_list, δs, L, n)
    create_plot_from_results(results, δs, plotdir)

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

    plotdir = dirname*"/parameter_sweep_g1"
    mkpath(plotdir)
    results = run_simulations(param_list, δs, L, n)
    create_plot_from_results(results, δs, plotdir)

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

    plotdir = dirname*"/parameter_sweep_g1_g2"
    mkpath(plotdir)
    results = run_simulations(param_list, δs, L, n)
    create_plot_from_results(results, δs, plotdir)

end

function main()

    # %% Manufactured solution
    L, m, r = (1, 1, 1)
    u_ex(x) = (x[1]^2 + x[2]^2  - 1)^2*100*cos(m*( 2π/L )*x[1])*cos(r*( 2π/L )*x[2])

    resultdir= "figures/biharmonic_CutCIP/"*string(Dates.now())
    println(resultdir)
    mkpath(resultdir)

    ns = [2^3, 2^4, 2^5, 2^6, 2^7]

    # @time convergence_analysis( ns=ns,  dirname=resultdir, u_ex=u_ex)
    @time translation_test(dirname=resultdir, u_ex=u_ex )

end

main()
