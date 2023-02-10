
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
    u_ex(x) = 1 - x[1]^2 - x[2]^2 -x[1]^3*x[2]
    ∇u_ex(x) = ∇(u_ex)(x)
    f(x) = - Δ(u_ex)(x)

    # First version
    # u_ex(x) = 1 - x[1]^2 - x[2]^2
    # ∇u_ex(x) = VectorValue(-2*x[1], -2*x[2])
    # f(x) = 4
    # ∇(::typeof(u_ex)) = ∇u_ex
    # ∇(u_ex) === ∇u_ex



    @with_kw struct Solution
        Ω
        Γ
        Ω_act
        Fg
        Λ

        bgmodel
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
        writevtk(sol.bgmodel,   vtkdirname*"/bgmodel")
        writevtk(sol.Ω,         vtkdirname*"/Omega")
        writevtk(sol.Ω_act,     vtkdirname*"/Omega_act")
        writevtk(sol.Λ,         vtkdirname*"/Lambda")
        writevtk(sol.Γ,         vtkdirname*"/Gamma")
        writevtk(sol.Fg,        vtkdirname*"/Fg")
    end


    function run(; order=order, n=n, vtkdirname=nothing )

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
        cutgeo_facets = cut_facets(bgmodel,geo)

        # Set up interpolation mesh and function spaces
        Ω_act = Triangulation(cutgeo, ACTIVE)

        # Construct function spaces
        V = TestFESpace(Ω_act, ReferenceFE(lagrangian, Float64, order), conformity=:L2)
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
        dΛ   = Measure(Λ, degree)
        dFg  = Measure(Fg, degree)

        # Set up normal vectors
        n_Γ = get_normal_vector(Γ)
        n_Λ = get_normal_vector(Λ)
        n_Fg = get_normal_vector(Fg)

        # Define weak form
        # Nitsche parameter
        # β = order*(order+1)
         β = 50


        # Ghost penalty parameter
        # γg = 0.1
        γg0 = β
        γg1 = 50
        γg2 = 0.1

        # Mesh size
        # h = (pmax - pmin)[1]/partition[1]
        h = L/n

        function mean_n(u,n)
            return 0.5*( u.plus⋅n.plus + u.minus⋅n.minus )
        end
        function jump_nn(u,n)
            return ( n.plus⋅ ∇∇(u).plus⋅ n.plus - n.minus ⋅ ∇∇(u).minus ⋅ n.minus )
            # return jump( n⋅ ∇∇(u)⋅ n)

        end
        # Define bilinear form
        a(u,v) =
            ∫( ∇(u)⋅∇(v) ) * dΩ  +
            # ∫( - jump(u)⋅mean_n(∇(v), n_Λ) - mean_n(∇(u), n_Λ)⋅jump(v) + (γd/h)*jump(u)⋅jump(v) ) * dΛ +  # version 1
            ∫( - jump(u ⋅ n_Λ)⋅mean(∇(v)) - mean(∇(u))⋅jump(v⋅n_Λ) + (β/h)*jump(u)⋅jump(v) ) * dΛ +      # version 2
            ∫( - u*(n_Γ⋅∇(v)) - (n_Γ⋅∇(u))*v + (β/h)*u*v ) * dΓ



        g(u,v) = ∫( (γg0/h)*jump(u)*jump(v)) * dFg +∫( (γg1*h)*jump(n_Fg⋅∇(u))*jump(n_Fg⋅∇(v)) ) * dFg

        if order == 2
            g(u,v) = ∫( (γg0/h)*jump(u)*jump(v)) * dFg +∫( (γg1*h)*jump(n_Fg⋅∇(u))*jump(n_Fg⋅∇(v)) ) * dFg +
                     ∫( (γg2*h^3)*jump_nn(u,n_Fg)*jump_nn(v,n_Fg) ) * dFg
        elseif order > 2
            println("Not supported order:", order)
        end

            # g(u,v)=  ∫( (γg1*h)*jump(n_Fg⋅∇(u))*jump(n_Fg⋅∇(v)) ) * dFg

        A(u,v) = a(u,v) + g(u,v)

        # Define linear form
        l(v) =
            ∫( f*v ) * dΩ +
            ∫( u_ex*(  - (n_Γ⋅∇(v))  + (β/h)*v  )) * dΓ

        # FE problem
        op = AffineFEOperator(A,l,U,V)
        uh = solve(op)

        e = u_ex - uh
        el2 = sqrt(sum( ∫(e*e)dΩ ))
        eh1 = sqrt(sum( ∫( e⊙e + ∇(e)⊙∇(e) )*dΩ ))
        eh_energy = sqrt(sum( ∫(∇(e)⊙∇(e) )*dΩ ))
        u_inter = interpolate(u_ex, V)

        sol = Solution(  bgmodel=bgmodel, Ω_act=Ω_act, Fg=Fg, Ω=Ω, Γ=Γ, Λ=Λ, h=h,
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

    resultdir= "figures/poisson_DGCutFEM/"*string(Dates.now())
    mkpath(resultdir)

    orders = [1,2]
    ns = [2^2, 2^3, 2^4, 2^5, 2^6]#, 2^7]
    dirname = resultdir
    makedir(dirname)
    convergence_analysis( orders=orders, ns=ns, dirname=dirname)
end

@time main()
