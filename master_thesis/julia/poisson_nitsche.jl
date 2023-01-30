include("results.jl")

module PoissonNitsche
    using Gridap
    using Parameters
    import GridapMakie
    import Makie
    import GLMakie
    using Test


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


    function run(; order=order, n=n, L=L,  m=m, r=r, use_quads=false)
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

        pmin = Point(0.,0.0)
        pmax = Point(L, L)
        partition = (n, n)

        if !use_quads
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
        reffe = ReferenceFE(lagrangian, Float64, order)

        V = TestFESpace(Ω, reffe, conformity=:H1)
        U = TrialFESpace(V)

        ## Define the weak form
        degree = 2*order
        dΩ = Measure(Ω, degree)
        dΓ= Measure(Γ, degree)

        n_Γ  = get_normal_vector(Γ)

        # Define mesh and stabilization parameters
        h = norm((pmax-pmin)./VectorValue(partition))

        # γ = 5.0*order*(order+1)  # Penalty parameter
        γ = 5  # Penalty parameter
        μ = γ/h

        a_Ω(u,v) =∫( ∇(v)⋅∇(u) )dΩ
        a_Γ(u,v) =∫( - ( ∇(u)⋅n_Γ )⊙v - u⊙( ∇(v)⋅n_Γ ) + μ*u⊙v )dΓ
        a(u,v) = a_Ω(u,v) + a_Γ(u,v)

        l_Ω(v) = ∫( v⊙f )dΩ
        l_Γ(v) = ∫( -(( ∇(v)⋅n_Γ )⊙g) + μ*(g⊙v) )dΓ
        l(v) = l_Ω(v) + l_Γ(v)


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







function convergence_analysis(; L, m, r, orders, ns, dirname, optimize)
    println("Run convergence",)

    for order in orders

        el2s = Float64[]
        eh1s = Float64[]
        ehs_energy = Float64[]
        println("Run convergence tests: order = "*string(order))

        for n in ns

            res = PoissonNitsche.run(order=order, n=n, L=L, m=m, r=r)

            if !(optimize)
                vtkdirname =dirname*"/order_"*string(order)*"_n_"*string(n)
                mkpath(vtkdirname)
                PoissonNitsche.generate_vtk(res=res, dirname=vtkdirname)
            end

            push!(el2s, res.el2)
            push!(eh1s, res.eh1)
            push!(ehs_energy, res.eh_energy)
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

    resultdir= "figures/poisson_nitsche/"*string(Dates.now())
    mkpath(resultdir)

    function run(;  L,m,r)
        orders=[2,3,4]
        ns = [2^2, 2^3, 2^4, 2^5]#, 2^6, 2^7]
        dirname = resultdir*"/L_"*string(round(L,digits=2))*"_m_"*string(m)*"_r_"*string(r);
        makedir(dirname)
        convergence_analysis( L=L, m=m, r=r, orders=orders, ns=ns, dirname=dirname, optimize=true)
    end
    run(L=1,m=1,r=1)
end

@time main()

