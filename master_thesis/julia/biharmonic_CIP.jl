
module BiharmonicCIP
    using Gridap
    using Parameters
    import GridapMakie
    import Makie
    import GLMakie
    using Test

    function man_sol(;L=1,m=1,r=1)
        u(x) = cos(m*( 2π/L )*x[1])*cos(r*( 2π/L )*x[2])
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


    function run(;order=order, n=n, L=L, m=m, r=r, use_quads=false)
        # Some parameters
        h = L/n
        γ = 1.5*order*( order+1)
        domain2D = (0, L, 0, L)
        partition2D = (n, n)
        if !use_quads
            model = CartesianDiscreteModel(domain2D,partition2D) |> simplexify
        else
            model = CartesianDiscreteModel(domain2D,partition2D)
        end

        # Manufactured solution
        function man_sol(;L=1,m=1,r=1)
            u(x) = cos(m*( 2π/L )*x[1])*cos(r*( 2π/L )*x[2])
        end
        u = man_sol(L=L, m=m, r=r)

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

        # manufactured solution
        # g(x) = ( ∇( Δ(u))⊙n_Γ)(x) #this does not compile since u is a ordinary function inner product with normal field vector
        α = 1
        f(x) = Δ(Δ(u))(x)+ α*u(x)
        # f(x) = ( 4 + α )*u(x)
        g(x) = 0

        function mean_nn(u,n)
            return 0.5*( n.plus⋅ ∇∇(u).plus⋅ n.plus + n.minus ⋅ ∇∇(u).minus ⋅ n.minus )
        end

        ⋅
        # Inner facets
        a(u,v) =( ∫( ∇∇(v)⊙∇∇(u) + α⋅(v⊙u) )dΩ
                 + ∫(-mean_nn(v,n_Λ)⊙jump(∇(u)⋅n_Λ) - mean_nn(u,n_Λ)⊙jump(∇(v)⋅n_Λ))dΛ
                 + ∫((γ/h)⋅jump(∇(u)⋅n_Λ)⊙jump(∇(v)⋅n_Λ))dΛ
                 + ∫(-( n_Γ ⋅ ∇∇(v)⋅ n_Γ )⊙∇(u)⋅n_Γ - ( n_Γ ⋅ ∇∇(u)⋅ n_Γ )⊙∇(v)⋅n_Γ)dΓ
                 + ∫((γ/h)⋅ ∇(u)⊙n_Γ⋅∇(v)⊙n_Γ )dΓ
                )

        l(v) = ∫( v ⋅ f )dΩ + ∫(- (g⋅v))dΓ

        op = AffineFEOperator(a, l, U, V)
        uh = solve(op)

        e = u - uh
        el2 = sqrt(sum( ∫(e*e)dΩ ))
        eh_energy = sqrt(sum( ∫( ∇∇(e)⊙∇∇(e) )*dΩ
                      + ( γ/h ) * ∫(jump(∇(e)⋅n_Λ) ⊙ jump(∇(e)⋅n_Λ))dΛ
                      + ( h/γ ) * ∫(mean_nn(e,n_Λ) ⊙ mean_nn(e,n_Λ))dΛ
                      + ( γ/h ) * ∫((∇(e)⋅n_Γ) ⊙ (∇(e)⋅n_Γ))dΓ
                      + ( h/γ ) * ∫(( n_Γ ⋅ ∇∇(e)⋅ n_Γ ) ⊙ ( n_Γ ⋅ ∇∇(e)⋅ n_Γ ))dΓ
                     ))

        eh1 = sqrt(sum( ∫( e⊙e + ∇(e)⊙∇(e) )*dΩ ))

        u_inter = interpolate(u, V)
        res = Results(  model=model, Ω=Ω, Γ=Γ, Λ=Λ, h=h,
                        u=u_inter, uh=uh, e=e, el2=el2, eh1=eh1, eh_energy=eh_energy)
        return res
    end

end # module



include("results.jl")




function convergence_analysis(; L, m, r, orders, ns, dirname, optimize)
    println("Run convergence",)

    for order in orders

        el2s = Float64[]
        eh1s = Float64[]
        ehs_energy = Float64[]
        println("Run convergence tests: order = "*string(order))

        for n in ns

            res = BiharmonicCIP.run(order=order, n=n, L=L, m=m, r=r)

            if !(optimize)
                vtkdirname =dirname*"/order_"*string(order)*"_n_"*string(n)
                mkpath(vtkdirname)
                BiharmonicCIP.generate_vtk(res=res, dirname=vtkdirname)
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

    resultdir= "figures/biharmonic_CIP/"*string(Dates.now())
    mkpath(resultdir)

    function run(;  L,m,r)
        orders=[2,3,4]
        ns = [2^2, 2^3, 2^4]#, 2^5]#, 2^6, 2^7]
        dirname = resultdir*"/L_"*string(round(L,digits=2))*"_m_"*string(m)*"_r_"*string(r);
        makedir(dirname)
        convergence_analysis( L=L, m=m, r=r, orders=orders, ns=ns, dirname=dirname, optimize=false)
    end

    run(L=1,m=1,r=1)
end

@time main()

