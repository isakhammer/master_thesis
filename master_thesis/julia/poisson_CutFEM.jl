
include("results.jl")
using Dates

module Solver
    using Gridap
    using GridapEmbedded
    using Parameters
    import Gridap: ∇


    # %% Manufactured solution
    # Provides a manufactured solution which is 0 on the unit circle
    # u_ex(x) = (x[1]^2 + x[2]^2  - 1)*sin(2π*x[1])*cos(2π*x[2])
    u_ex(x) = 1 - x[1]^2 - x[2]^2
    f(x) = 4
    ∇u_ex(x) = VectorValue(-2*x[1], -2*x[2])

    ∇(::typeof(u_ex)) = ∇u_ex
    ∇(u_ex) === ∇u_ex

    # function man_sol(;L=1,m=1,r=1)
    #     u(x) = 100*cos(m*( 2π/L )*x[1])*cos(r*( 2π/L )*x[2])
    # end


    @with_kw struct Solution
        Ω
        Γ
        Ω_act
        Fg

        bgmodel
        h::Real

        u
        uh
        e
        el2
        eh1
        eh_energy
    end

    function generate_vtk(; res::Solution, dirname::String)
        println("Generating vtk's in ", dirname)
        if (isdir(dirname))
            rm(dirname, recursive=true)
        end
        mkdir(dirname)

        # Write out models and computational domains for inspection
        writevtk(res.bgmodel, dirname*"/bgmodel")
        writevtk(res.Ω, dirname*"/Omega")
        writevtk(res.Ω_act, dirname*"/Omega_act")
        writevtk(res.Γ, dirname*"/Gamma")
        writevtk(res.Fg, dirname*"/Fg")

    end


    function run(; order=order, n=n, dirname=nothing )

        # Background model
        L = 1.11
        domain = (-L, L, -L, L)
        pmin = Point(-L, -L)
        pmax = Point(L, L)
        partition = (n,n)
        bgmodel = CartesianDiscreteModel(pmin, pmax, partition)

        # Implicit geometry
        # TODO: Define own level set function via AnalyticalGeometry
        R  = 1.0
        geo = disk(R)

        # Cut the background model
        cutgeo = cut(bgmodel, geo)

        # Set up interpolation mesh and function spaces
        Ω_act = Triangulation(cutgeo, ACTIVE)

        # Construct function spaces
        V = TestFESpace(Ω_act, ReferenceFE(lagrangian, Float64, order),conformity=:H1)
        U = TrialFESpace(V)

        # Set up integration meshes, measures and normals
        Ω = Triangulation(cutgeo, PHYSICAL)
        Γ = EmbeddedBoundary(cutgeo)
        Fg = GhostSkeleton(cutgeo)

        # Set up integration measures
        degree = 2*order
        dΩ   = Measure(Ω, degree)
        dΓ   = Measure(Γ, degree)
        dFg  = Measure(Fg, degree)

        # Set up normal vectors
        n_Γ = get_normal_vector(Γ)
        n_Fg = get_normal_vector(Fg)


        # Define weak form
        # Nitsche parameter
        γd = order(order+1)

        # Ghost penalty parameter
        γg = 0.1

        # Mesh size
        # h = (pmax - pmin)[1]/partition[1]
        h = L/n

        # Define bilinear form
        a(u,v) =
            ∫( ∇(u)⋅∇(v) ) * dΩ  +
            ∫( (γd/h)*u*v  - u*(n_Γ⋅∇(v)) - (n_Γ⋅∇(u))*v ) * dΓ +
            ∫( (γg*h)*jump(n_Fg⋅∇(u))*jump(n_Fg⋅∇(v)) ) * dFg

        # Define linear form
        l(v) =
            ∫( f*v ) * dΩ +
            ∫( u_ex*( (γg/h)*v - (n_Γ⋅∇(v)) )  ) * dΓ

        # FE problem
        op = AffineFEOperator(a,l,U,V)
        uh = solve(op)


        e = u_ex - uh

        el2 = sqrt(sum( ∫(e*e)dΩ ))
        eh1 = sqrt(sum( ∫( e⊙e + ∇(e)⊙∇(e) )*dΩ ))
        eh_energy = sqrt(sum( ∫(∇(e)⊙∇(e) )*dΩ ))
        u_inter = interpolate(u_ex, V)

        sol = Solution(  bgmodel=bgmodel, Ω_act=Ω_act, Fg=Fg, Ω=Ω, Γ=Γ,  h=h,
                        u=u_inter, uh=uh, e=e, el2=el2, eh1=eh1, eh_energy=eh_energy)

        if ( dirname!=nothing)
            generate_vtk(sol, dirname)
        end

        return sol
    end

end # module



function convergence_analysis(; orders, ns, dirname, optimize=true)
    println("Run convergence",)

    for order in orders

        el2s = Float64[]
        eh1s = Float64[]
        ehs_energy = Float64[]
        println("Run convergence tests: order = "*string(order))

        for n in ns

            sol = nothing
            if !(optimize)
                vtkdirname =dirname*"/order_"*string(order)*"_n_"*string(n)
                mkpath(vtkdirname)
                sol = Solver.run(order=order, n=n, dirname=vtkdirname)
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

    resultdir= "figures/poisson_CutFEM/"*string(Dates.now())
    mkpath(resultdir)

    function run()
        orders=[1,2,3,4]
        ns = [2^2, 2^3, 2^4, 2^5, 2^6, 2^7]
        dirname = resultdir
        makedir(dirname)
        convergence_analysis( orders=orders, ns=ns, dirname=dirname)
    end

    run()
end

@time main()

