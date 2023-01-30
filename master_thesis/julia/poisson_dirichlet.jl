
module PoissonDirichlet
    using Gridap
    using Parameters
    import GridapMakie
    import Makie
    import GLMakie
    using Test

    @with_kw struct SolverSettings

        order::Int      # Order on elements

        # Domain Specific
        L::Real         # Square length
        n::Int          # Number of partitions
        use_quads::Bool = true

        # Manufactured solution parameters
        m::Int
        r::Int
    end

    @with_kw struct Results
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
    end

    function generate_vtk(; res::Results, dirname::String)
        println("Generating vtk's in ", dirname)
        if (isdir(dirname))
            rm(dirname, recursive=true)
        end
        mkdir(dirname)

        writevtk(res.model, dirname*"/model")
        writevtk(res.Λ, dirname*"/skeleton")
        writevtk(res.Γ, dirname*"/boundary")
        writevtk(res.Λ, dirname*"/jumps",cellfields=["jump_u"=>jump(res.uh)])
        writevtk(res.Ω, dirname*"/omega",cellfields=["uh"=>res.uh])
        writevtk(res.Ω, dirname*"/error",cellfields=["e"=>res.e])
        writevtk(res.Ω, dirname*"/manufatured",cellfields=["u"=>res.u])


        fig = Makie.plot(res.Λ)
        Makie.wireframe!(res.Λ, color=:black, linewidth=2)
        Makie.wireframe!(res.Γ, color=:black, linewidth=2)
        Makie.save(dirname*"/grid.png", fig)

        # (Isak): Doesnt work :( Please fix
        # fig, _ , plt = Makie.plot(res.Ω, res.uh)
        # Makie.Colorbar(fig[1,2], plt)
        # Makie.save(dirname*"/man_sol.png", fig)
    end


    function run(set::SolverSettings)
        # Some parameters
        # h = set.L/set.n
        # γ = 1.5*set.order*( set.order+1)
        # domain2D = (0, set.L, 0, set.L)
        # partition2D = (set.n, set.n)
        # if !set.use_quads
        #     model = CartesianDiscreteModel(domain2D,partition2D) |> simplexify
        # else
        #     model = CartesianDiscreteModel(domain2D,partition2D)
        # end

        ##
        n = set.n

        pmin = Point(0.,0.0)
        pmax = Point(2π, 2π)
        partition = (n, n)

        if !set.use_quads
            model = CartesianDiscreteModel(pmin, pmax, partition) |> simplexify
        else
            model = CartesianDiscreteModel(pmin, pmax, partition)
        end
        h = norm((pmax-pmin)./VectorValue(partition))
        ##

        # u is the manufactured solution
        # -Δu = f in Ω, and u = g on Γ
        # u = man_sol(L=set.L, m=set.m, r=set.r)
        u(x) = 100*cos(x[1])*cos(x[2])
        f(x) = 2*u(x)
        g(x) = u(x)

        # Define triangulation
        Ω = Triangulation(model)
        Γ = BoundaryTriangulation(model)
        Λ = SkeletonTriangulation(model)

        ## Function spaces
        # order = 2
        # reffe = ReferenceFE(lagrangian, Float64, order)
        reffe = ReferenceFE(lagrangian, Float64, set.order)

        V = TestFESpace(Ω, reffe, conformity=:H1, dirichlet_tags="boundary")
        U = TrialFESpace(V,g)

        ## Define the weak form
        degree = 2*set.order
        dΩ = Measure(Ω, degree)

        a_Ω(u,v) =∫( ∇(v)⋅∇(u) )dΩ
        a(u,v) = a_Ω(u,v)

        l_Ω(v) = ∫( v⊙f )dΩ
        l(v) = l_Ω(v)

        op = AffineFEOperator(a, l, U, V)
        uh = solve(op)

        e = u - uh
        el2 = sqrt(sum( ∫(e*e)dΩ ))
        eh_energy = sqrt(sum( ∫(∇(e)⊙∇(e) )*dΩ ))
        eh1 = sqrt(sum( ∫( e⊙e + ∇(e)⊙∇(e) )*dΩ ))

        u_inter = interpolate(u, V)
        res = Results(  model=model, Ω=Ω, Γ=Γ, Λ=Λ, h=h,
                        u=u_inter, uh=uh, e=e, el2=el2, eh1=eh1, eh_energy=eh_energy)
        return res
    end

end # module




