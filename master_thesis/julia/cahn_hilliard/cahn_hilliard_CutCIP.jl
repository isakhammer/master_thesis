
using Test
using Gridap
using Plots
using Dates
using LaTeXStrings

module Solver
    using Gridap
    using LinearAlgebra
    using PROPACK
    using GridapEmbedded
    using GridapGmsh
    using Parameters
    import Gridap: ∇

    # α(x) = x[1]^2 + x[2]^2
    α = 1

    function man_sol(u_ex)
        α = 1
        f(t) = x ->  Δ(u_ex(t))(x)
        ∇u_ex(t) = x ->  ∇(u_ex(t))(x)
        ∇Δu_ex(t) = x ->  ∇(Δ(u_ex(t)))(x)
        return u_ex, f, ∇u_ex, ∇Δu_ex
    end

    @with_kw struct Config
        exact_sol
        circle
    end

    function run(;n, dt, solver_config,  vtkdirname=nothing)
        u_ex, f, ∇u_ex, ∇Δu_ex = solver_config.exact_sol
        order=2

        # Background model
        L = 1.11
        domain = (-L, L, -L, L)
        pmin = Point(-L, -L)
        pmax = Point(L, L)
        partition = (n,n)
        bgmodel = CartesianDiscreteModel(pmin, pmax, partition)

        # Implicit geometry
        R  = 1.0
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

        # Define weak form
        γ = 10

        # Ghost penalty parameter
        # γg0 = γ
        γg1 = 10/2
        γg2 = 0.1

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


        # Define weak form
        γ = 10

        # Ghost penalty parameter
        # γg0 = γ
        γg1 = 10/2
        γg2 = 0.1

        # Inner facets
        a(t,u,v) =( ∫( ∇∇(v)⊙∇∇(u) + α⋅(v⊙u) )dΩ
                 + ∫(-mean_nn(v,n_Λ)⊙jump(∇(u)⋅n_Λ) - mean_nn(u,n_Λ)⊙jump(∇(v)⋅n_Λ))dΛ
                 + ∫(-( n_Γ ⋅ ∇∇(v)⋅ n_Γ )⊙∇(u)⋅n_Γ - ( n_Γ ⋅ ∇∇(u)⋅ n_Γ )⊙∇(v)⋅n_Γ)dΓ
                 + ∫((γ/h)⋅jump(∇(u)⋅n_Λ)⊙jump(∇(v)⋅n_Λ))dΛ + ∫((γ/h)⋅ ∇(u)⊙n_Γ⋅∇(v)⊙n_Γ )dΓ
                )

        # Define linear form
        g_1(t) = ∇u_ex(t)⋅n_Γ
        g_2(t) = ∇Δu_ex(t)⋅n_Γ
        l(t,v) = (∫( f(t)*v ) * dΩ
                  +  ∫(-(g_2(t)⋅v))dΓ
                  + ∫(g_1(t)⊙(-(n_Γ⋅∇∇(v)⋅n_Γ) + (γ/h)*∇(v)⋅n_Γ)) * dΓ
               )

        # # Define of ghost penalties
        if order != 2
            println("Not supported order:", order)
        end

        g(t,u,v) = h^(-2)*( ∫( (γg1*h)*jump(n_Fg⋅∇(u))*jump(n_Fg⋅∇(v)) ) * dFg +
                         ∫( (γg2*h^3)*jump_nn(u,n_Fg)*jump_nn(v,n_Fg) ) * dFg)

        A(t,u,v) = a(t,u,v) + g(t,u,v)

        # # bilinear form
        m(t, u, v) = ∫( α* u⋅v )dΩ

    end

end # Solver

function convergence_analysis(; n, dts, dirname, solver_config, write_vtks=true)
    println("Run convergence",)

    el2s = Float64[]
    eh1s = Float64[]
    ehs_energy = Float64[]
    cond_numbers = Float64[]
    ndofs = Float64[]
    println("Run convergence tests: n = "*string(n))

    for dt in dts
        sol = Solver.run(n=n, dt=dt, solver_config=solver_config, vtkdirname=dirname)
        #     if (write_vtks)
        #     else
        #         sol = Solver.run(order=order, solver_config, n=n)
        #     end

        #     push!(el2s, sol.el2)
        #     push!(eh1s, sol.eh1)
        #     push!(ehs_energy, sol.eh_energy)
        #     push!(cond_numbers, sol.cond_number)
        #     push!(ndofs, sol.ndof)
        # end
        # generate_figures(ns=ns, el2s=el2s, eh1s=eh1s, ehs_energy=ehs_energy,
        #                  cond_numbers=cond_numbers, order=order, ndofs=ndofs, dirname=dirname)
    end
end

function main()

    dirname= "figures/CIP_cahn_hilliard/example"*string(Dates.now())
    println(dirname)
    mkpath(dirname)

    u_ex(x,t::Real) = t*cos(x[1])*cos(x[2])
    u_ex(t::Real) = x -> u_ex(x,t)
    exact_sol = Solver.man_sol(u_ex)
    circle = true
    solver_config = Solver.Config(exact_sol, circle)

    dts = [2^-4]
    @time convergence_analysis( n=2^6, dts=dts, solver_config=solver_config, dirname=dirname)

end

main()
