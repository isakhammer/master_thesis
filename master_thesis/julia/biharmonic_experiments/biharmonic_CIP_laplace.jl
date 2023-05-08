using Dates
import Plots
Plots.pyplot()

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

    function compute_condition_number(A, ndof; method="p-inf")

        if method == "svd"
            tolin = 10^-8
            s_max, _ = tsvdvals(A, k = 1, tolin = tolin)
            s_min, _ = tsvdvals_irl(A, k = 2, tolin = tolin)
            # Take the smallest singular value which isn't 0
            s_min = isapprox(s_min[2], 0.0, atol=10^-10) ? s_min[1] : s_min[2]
            return s_max[1] / s_min
        elseif method == "p-inf"
            # https://github.com/mfasi/julia/blob/4ceb4ea9ee46ea92d406cfced308451a112d16f9/base/sparse/linalg.jl#L505
            return ( 1/sqrt(ndof) )*cond(A,Inf)
        else method == "p-1"
            return cond(A,p=1)
        end

        # catch e
        #     # Use more expansive full svd if PROPACK throw "Invariant Subspace" error
        #     svd_values = svdvals(Matrix(A))
        #     return svd_values[1]/svd_values[end]
        # end
    end


    function man_sol(u_ex)
        ∇u_ex(x) = ∇(u_ex)(x)
        ∇Δu_ex(x) = ∇(Δ(u_ex))(x)
        f(x) = Δ(Δ(u_ex))(x)+ α(x)⋅u_ex(x)
        return u_ex, f, ∇u_ex, ∇Δu_ex
    end

    @with_kw struct Config
        exact_sol
        circle
    end

    @with_kw struct Solution
        Ω
        Γ
        Λ

        model
        h::Real

        u
        uh
        e
        el2
        eh1
        eh_energy
        cond_number
        ndof
    end

    function generate_vtk(; sol::Solution, dirname::String)

        # Write out models and computational domains for inspection
        writevtk(sol.model,   dirname*"/model")
        writevtk(sol.Ω,         dirname*"/Omega")
        # writevtk(sol.Ω_act,   dirname*"/Omega_act")
        writevtk(sol.Λ,         dirname*"/Lambda")
        writevtk(sol.Γ,         dirname*"/Gamma")
        # writevtk(sol.Fg,        dirname*"/Fg")
        writevtk(sol.Λ,         dirname*"/jumps",        cellfields=["jump_u"=>jump(sol.uh)])
        writevtk(sol.Ω,         dirname*"/omega",        cellfields=["uh"=>sol.uh])
        writevtk(sol.Ω,         dirname*"/error",        cellfields=["e"=>sol.e])
        writevtk(sol.Ω,         dirname*"/manufatured",  cellfields=["u"=>sol.u])
    end


    function run(;order, n, solver_config, vtkdirname=nothing)

        u_ex, f, ∇u_ex, ∇Δu_ex = solver_config.exact_sol
        if solver_config.circle==true
            refinement = Integer(log2(n)) - 1
            if ( refinement < 1 || refinement > 6 )
                println("Refinement out of range.")
            end
            mesh_file = "circle/circle_"*string(refinement)*".msh"
            model = GmshDiscreteModel(mesh_file)
            h = 2/n

        else
            # Some parameters
            L=1
            h = L/n
            domain2D = (0, L, 0, L)
            partition2D = (n, n)
            model = CartesianDiscreteModel(domain2D,partition2D)
            # model = CartesianDiscreteModel(domain2D,partition2D) |> simplexify
        end

        γ = 1.5*order*( order+1)

        # Spaces
        V = TestFESpace(model, ReferenceFE(lagrangian,Float64, order), conformity=:H1)
        U = TrialFESpace(V)
        Ω = Triangulation(model)
        Λ = SkeletonTriangulation(model)
        Γ = BoundaryTriangulation(model)

        degree = 2*order
        dΩ = Measure(Ω,degree)
        dΓ = Measure(Γ,degree)
        dΛ = Measure(Λ,degree)

        n_Λ = get_normal_vector(Λ)
        n_Γ = get_normal_vector(Γ)

        function mean_nn(u,n)
            return 0.5*( n.plus⋅ ∇∇(u).plus⋅ n.plus + n.minus ⋅ ∇∇(u).minus ⋅ n.minus )
        end

        # Inner facets
        a(u,v) =( ∫( Δ(v)⊙Δ(u) + α⋅(v⊙u) )dΩ
                 + ∫(-mean(Δ(v))⊙jump(∇(u)⋅n_Λ) - mean(Δ(u))⊙jump(∇(v)⋅n_Λ) + (γ/h)⋅jump(∇(u)⋅n_Λ)⊙jump(∇(v)⋅n_Λ))dΛ
                 + ∫(-Δ(v)⊙∇(u)⋅n_Γ - Δ(u)⊙∇(v)⋅n_Γ + (γ/h)⋅ ∇(u)⊙n_Γ⋅∇(v)⊙n_Γ )dΓ
                )

        g_1 = ∇u_ex⋅n_Γ
        g_2 = ∇Δu_ex⋅n_Γ

        l(v) = (∫( f*v ) * dΩ
                +  ∫(-(g_2⋅v))dΓ
                + ∫(g_1⊙(-Δ(v) + (γ/h)*∇(v)⋅n_Γ)) * dΓ
               )

        op = AffineFEOperator(a, l, U, V)
        uh = solve(op)
        A =  get_matrix(op)
        ndof = size(A)[1]
        h = 1/sqrt(ndof)
        cond_number = compute_condition_number(A, ndof)

        e = u_ex - uh
        el2 = sqrt(sum( ∫(e*e)dΩ ))

        # TODO: Add α into ∫(e⊙e)*dΩ
        eh_energy = sqrt(sum( ∫(e⊙e)*dΩ + ∫( ∇∇(e)⊙∇∇(e) )*dΩ
                      + ( γ/h ) * ∫(jump(∇(e)⋅n_Λ) ⊙ jump(∇(e)⋅n_Λ))dΛ
                      + ( h/γ ) * ∫(mean_nn(e,n_Λ) ⊙ mean_nn(e,n_Λ))dΛ
                      + ( γ/h ) * ∫((∇(e)⋅n_Γ) ⊙ (∇(e)⋅n_Γ))dΓ
                      + ( h/γ ) * ∫(( n_Γ ⋅ ∇∇(e)⋅ n_Γ ) ⊙ ( n_Γ ⋅ ∇∇(e)⋅ n_Γ ))dΓ
                     ))

        eh1 = sqrt(sum( ∫( e⊙e + ∇(e)⊙∇(e) )*dΩ ))

        u_inter = interpolate(u_ex, V)


        sol = Solution(  model=model, Ω=Ω, Γ=Γ, Λ=Λ, h=h,
                        u=u_inter, uh=uh, e=e, el2=el2, eh1=eh1, eh_energy=eh_energy,
                        cond_number=cond_number, ndof=ndof)

        if ( vtkdirname!=nothing)
            dirname =vtkdirname*"/order_"*string(order)*"_n_"*string(n)
            mkpath(dirname)
            println("Generating vtk's in ", dirname)

            generate_vtk(sol=sol, dirname=dirname)
        end

        return sol
    end

end # module


function generate_figures(;ns, el2s, eh1s, ehs_energy, cond_numbers, ndofs, order::Integer, dirname::String)
    filename = dirname*"/conv_order_"*string(order)

    hs = 1 .// ns
    # hs_str =  latexify.(hs)

    compute_eoc(hs, errs) = log.(errs[1:end-1]./errs[2:end])./log.(hs[1:end-1]./hs[2:end])
    eoc_l2 = compute_eoc(hs, el2s)
    eoc_eh1 = compute_eoc(hs,eh1s)
    eoc_eh_energy = compute_eoc(hs,ehs_energy)

    eoc_l2 =  [nothing; eoc_l2]
    eoc_eh1 =  [nothing; eoc_eh1]
    eoc_eh_energy =  [nothing; eoc_eh_energy]
    data = hcat(ns, el2s,  eoc_l2, eh1s, eoc_eh1, ehs_energy, eoc_eh_energy, cond_numbers)
    header = [L"$n$", L"$\Vert e \Vert_{L^2}$", "EOC", L"$ \Vert e \Vert_{H^1}$", "EOC", L"$\Vert e \Vert_{ a_h,* }$", "EOC", "Cond number"]

    open(filename*".tex", "w") do io
        pretty_table(io, data, header=header, backend=Val(:latex ), formatters = ( ft_printf("%.3E"), ft_nonothing ))
    end

    minimal_header = ["n", "L2", "EOC", "H1", "EOC", "a_h", "EOC", "cond", "const", "ndofs"]
    data = hcat(ns, el2s,  eoc_l2, eh1s, eoc_eh1, ehs_energy, eoc_eh_energy, cond_numbers, cond_numbers.*hs.^4, ndofs)
    pretty_table(data, header=minimal_header, formatters = ( ft_printf("%.0f",[1,10]), ft_printf("%.2f",[3,5,7]), ft_printf("%.1E",[2,4,6,8,9]), ft_nonothing ))

    # Plots.plot(hs, (el2s, eh1s, ehs_energy), label=(L"\\Vert e \\Vert_{L^2}", L"\\Vert e \\Vert_{ H^1 }^{  } ", L"\\Vert e  \\Vert_{ a_h,* }"))
    Plots.plot(hs, [el2s, eh1s, ehs_energy], label=[L" | e |_{L^2}", L"| e |_{ H^1 }^{  } ", L"| e  |_{ a_h,* }"])
    Plots.plot(hs, el2s, label=L" | e |_{L^2}")
    Plots.plot!(hs, eh1s, label=L"| e |_{ H^1 }^{  } ")
    Plots.plot!(hs, ehs_energy, label=L"| e  |_{ a_h,* }")

    Plots.plot!(xscale=:log2, yscale=:log2, minorgrid=true)
    Plots.scatter!(hs, [el2s, eh1s, ehs_energy],primary=false)
    Plots.xlabel!("h")
    Plots.savefig(filename*".png")
end

function convergence_analysis(; orders, ns, dirname, solver_config, write_vtks=true)
    println("Run convergence",)

    for order in orders

        el2s = Float64[]
        eh1s = Float64[]
        ehs_energy = Float64[]
        cond_numbers = Float64[]
        ndofs = Float64[]
        println("Run convergence tests: order = "*string(order))

        for n in ns

            if (write_vtks)
                sol = Solver.run(order=order, n=n, solver_config=solver_config, vtkdirname=dirname)
            else
                sol = Solver.run(order=order, solver_config, n=n)
            end

            push!(el2s, sol.el2)
            push!(eh1s, sol.eh1)
            push!(ehs_energy, sol.eh_energy)
            push!(cond_numbers, sol.cond_number)
            push!(ndofs, sol.ndof)
        end
        generate_figures(ns=ns, el2s=el2s, eh1s=eh1s, ehs_energy=ehs_energy,
                         cond_numbers=cond_numbers, order=order, ndofs=ndofs, dirname=dirname)
    end
end


function main()

    # %% Manufactured solution
    L, m, r = (1, 1, 1)
    # u_ex(x) = (x[1]^2 + x[2]^2  - 1)^2*sin(2π*x[1])*cos(2π*x[2])
    # u_ex(x) = 100*sin(m*( 2π/L )*x[1])*cos(r*( 2π/L )*x[2])
    u_ex(x) = 100*cos(m*( 2π/L )*x[1])*cos(r*( 2π/L )*x[2])
    exact_sol = Solver.man_sol(u_ex)
    circle = false
    solver_config = Solver.Config(exact_sol, circle)

    resultdir= "figures/biharmonic_CIP_nitsche/"*string(Dates.now())
    println(resultdir)
    mkpath(resultdir)

    orders = [2]
    # orders = [2,3,4]
    # ns = [2^2, 2^3, 2^4, 2^5]
    ns = [2^2, 2^3, 2^4, 2^5, 2^6, 2^7]
    dirname = resultdir
    mkpath(dirname)
    @time convergence_analysis( orders=orders, ns=ns, solver_config=solver_config, dirname=dirname)
end


main()

