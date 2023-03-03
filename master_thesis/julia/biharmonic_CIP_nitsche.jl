
include("results.jl")
using Dates

module Solver
    using Gridap
    using GridapEmbedded
    using Parameters
    import Gridap: ∇


    # %% Manufactured solution
    # Provides a manufactured solution which is 0 on the unit circle

    L, m, r = (1, 1, 1)
    # u_ex(x) = (x[1]^2 + x[2]^2  - 1)*sin(2π*x[1])*cos(2π*x[2])
    u_ex(x) = 100*cos(m*( 2π/L )*x[1])*cos(r*( 2π/L )*x[2])

    # u_ex(x) = 1 - x[1]^2 - x[2]^2 -x[1]^3*x[2]
    ∇u_ex(x) = ∇(u_ex)(x)
    ∇Δu_ex(x) = ∇(Δ(u_ex))(x)

    α = 1
    f(x) = Δ(Δ(u_ex))(x)+ α*u_ex(x)

    @with_kw struct Solution
        Ω
        Γ
        # Ω_act
        # Fg
        Λ

        model
        h::Real

        u
        uh
        e
        el2
        eh1
        eh_energy
    end

    function generate_vtk(; sol::Solution, vtkdirname::String)
        println("Generating vtk's in ", vtkdirname)
        mkpath(vtkdirname)

        # Write out models and computational domains for inspection
        writevtk(sol.model,   vtkdirname*"/model")
        writevtk(sol.Ω,         vtkdirname*"/Omega")
        # writevtk(sol.Ω_act,     vtkdirname*"/Omega_act")
        writevtk(sol.Λ,         vtkdirname*"/Lambda")
        writevtk(sol.Γ,         vtkdirname*"/Gamma")
        # writevtk(sol.Fg,        vtkdirname*"/Fg")
        writevtk(sol.Λ,         vtkdirname*"/jumps",        cellfields=["jump_u"=>jump(sol.uh)])
        writevtk(sol.Ω,         vtkdirname*"/omega",        cellfields=["uh"=>sol.uh])
        writevtk(sol.Ω,         vtkdirname*"/error",        cellfields=["e"=>sol.e])
        writevtk(sol.Ω,         vtkdirname*"/manufatured",  cellfields=["u"=>sol.u])
    end


    function run(;order=order, n=n, vtkdirname=nothing)
        # Some parameters
        h = L/n
        γ = 1.5*order*( order+1)
        domain2D = (0, L, 0, L)
        partition2D = (n, n)
        use_quads=false
        if !use_quads
            model = CartesianDiscreteModel(domain2D,partition2D) |> simplexify
        else
            model = CartesianDiscreteModel(domain2D,partition2D)
        end

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
        a(u,v) =( ∫( ∇∇(v)⊙∇∇(u) + α⋅(v⊙u) )dΩ
                 + ∫(-mean_nn(v,n_Λ)⊙jump(∇(u)⋅n_Λ) - mean_nn(u,n_Λ)⊙jump(∇(v)⋅n_Λ))dΛ
                 + ∫(-( n_Γ ⋅ ∇∇(v)⋅ n_Γ )⊙∇(u)⋅n_Γ - ( n_Γ ⋅ ∇∇(u)⋅ n_Γ )⊙∇(v)⋅n_Γ)dΓ
                 + ∫((γ/h)⋅jump(∇(u)⋅n_Λ)⊙jump(∇(v)⋅n_Λ))dΛ + ∫((γ/h)⋅ ∇(u)⊙n_Γ⋅∇(v)⊙n_Γ )dΓ
                )

        g_1 = ∇u_ex⋅n_Γ
        g_2 = ∇Δu_ex⋅n_Γ

        l(v) = (∫( f*v ) * dΩ
                +  ∫(-(g_2⋅v))dΓ
                + ∫(g_1⊙(-(n_Γ⋅∇∇(v)⋅n_Γ) + (γ/h)*∇(v)⋅n_Γ)) * dΓ
               )

        op = AffineFEOperator(a, l, U, V)
        uh = solve(op)

        e = u_ex - uh
        el2 = sqrt(sum( ∫(e*e)dΩ ))
        eh_energy = sqrt(sum( ∫( ∇∇(e)⊙∇∇(e) )*dΩ
                      + ( γ/h ) * ∫(jump(∇(e)⋅n_Λ) ⊙ jump(∇(e)⋅n_Λ))dΛ
                      + ( h/γ ) * ∫(mean_nn(e,n_Λ) ⊙ mean_nn(e,n_Λ))dΛ
                      + ( γ/h ) * ∫((∇(e)⋅n_Γ) ⊙ (∇(e)⋅n_Γ))dΓ
                      + ( h/γ ) * ∫(( n_Γ ⋅ ∇∇(e)⋅ n_Γ ) ⊙ ( n_Γ ⋅ ∇∇(e)⋅ n_Γ ))dΓ
                     ))

        eh1 = sqrt(sum( ∫( e⊙e + ∇(e)⊙∇(e) )*dΩ ))

        u_inter = interpolate(u_ex, V)
        sol = Solution(  model=model, Ω=Ω, Γ=Γ, Λ=Λ, h=h,
                        u=u_inter, uh=uh, e=e, el2=el2, eh1=eh1, eh_energy=eh_energy)

        if ( vtkdirname!=nothing)
            generate_vtk(sol=sol, vtkdirname=vtkdirname)
        end

        return sol
    end

end # module


function convergence_analysis(; orders, ns, dirname, write_vtks=true)
    println("Run convergence",)

    for order in orders

        el2s = Float64[]
        eh1s = Float64[]
        ehs_energy = Float64[]
        println("Run convergence tests: order = "*string(order))

        for n in ns

            if (write_vtks)
                vtkdirname =dirname*"/order_"*string(order)*"_n_"*string(n)
                mkpath(vtkdirname)
                sol = Solver.run(order=order, n=n, vtkdirname=vtkdirname)
            else
                sol = Solver.run(order=order, n=n)
            end

            push!(el2s, sol.el2)
            push!(eh1s, sol.eh1)
            push!(ehs_energy, sol.eh_energy)
        end
        Results.generate_figures(ns=ns, el2s=el2s, eh1s=eh1s, ehs_energy=ehs_energy, order=order, dirname=dirname)
    end
end


function main()
    function makedir(dirname)
        if (isdir(dirname))
            rm(dirname, recursive=true)
        end
        mkdir(dirname)
    end

    resultdir= "figures/biharmonic_CIP_nitsche/"*string(Dates.now())
    println(resultdir)
    mkpath(resultdir)

    orders = [2,3,4]
    ns = [2^2, 2^3, 2^4, 2^5, 2^6]#, 2^7]
    dirname = resultdir
    makedir(dirname)
    convergence_analysis( orders=orders, ns=ns, dirname=dirname)
end

@time main()
