
module BiharmonicEquation
    using Gridap
    using Parameters
    using GridapMakie
    using GLMakie
    using Test


    @with_kw struct Results
        Ω
        Γ
        Λ

        model
        h::Real
        γ::Real
        order::Int
        degree::Int

        u_inter
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

        # writevtk(res.model, dirname*"/model")
        # writevtk(res.Λ, dirname*"/skeleton")
        # writevtk(res.Γ, dirname*"/boundary")
        # writevtk(res.Λ, dirname*"/jumps",cellfields=["jump_u"=>jump(res.uh)])
        # writevtk(res.Ω, dirname*"/omega",cellfields=["uh"=>res.uh])
        # writevtk(res.Ω, dirname*"/error",cellfields=["e"=>res.e])
        # writevtk(res.Ω, dirname*"/manufatured",cellfields=["u"=>res.u])


        fig = plot(res.Λ)
        wireframe!(res.Λ, color=:black, linewidth=2)
        wireframe!(res.Γ, color=:black, linewidth=2)
        save(dirname*"/grid.png", fig)

        fig, _ , plt = plot(res.Ω, res.u_inter)
        Colorbar(fig[1,2], plt)
        save(dirname*"/man_sol.png", fig)

    end


    function run_CP_method(;n, L, γ, order,  u::Function, simplex=true)
        h = L/n
        domain2D = (0, L, 0, L)
        partition2D = (n,n)
        model = CartesianDiscreteModel(domain2D,partition2D)

        # Spaces
        V = TestFESpace(model, ReferenceFE(lagrangian,Float64,order), conformity=:H1)
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

        res = Results( model=model, Ω=Ω, Γ=Γ, Λ=Λ,
                      h=h, γ=γ, order=order, degree=degree,
                      u_inter=u, uh=uh, e=e, el2=el2, eh1=eh1, eh_energy=eh_energy)

        return res
    end


    function man_sol(;L=1,m=1,r=1)
        u(x) = cos(m*( 2π/L )*x[1])*cos(r*( 2π/L )*x[2])
    end

end # module






