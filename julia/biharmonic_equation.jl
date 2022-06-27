
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
        eh
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


    function run_CP_method(;n, L, γ, order,  u::Function = man_sol(), simplex=true)
        h = L / n
        γ = 5*order*(order+1)

        domain = (0,L,0,L)
        partition = (n,n)

        if simplex
            model = CartesianDiscreteModel(domain,partition) |> simplexify
        else
            model = CartesianDiscreteModel(domain,partition)
        end

        # FE space
        V = TestFESpace(model,ReferenceFE(lagrangian,Float64,order),conformity=:H1)
        U = TrialFESpace(V)

        # Triangulation
        Ω = Triangulation(model)
        Γ = BoundaryTriangulation(model)
        Λ = SkeletonTriangulation(model)
        degree = 2*order
        dΩ = Measure(Ω,degree)
        dΓ = Measure(Γ,degree)
        dΛ = Measure(Λ,degree)
        n_Γ = get_normal_vector(Γ)
        n_Λ = get_normal_vector(Λ)

        # Analytical manufactured solution
        α = 1
        f(x) = Δ(Δ(u))(x)+ α*u(x)
        g(x) = Δ(u)(x)

        # PROBLEM 1
        # Statement: laplacian Δ(u) makes no sense since we need the hessian
        # we should use ∇∇(u) instead:
        # Definition of  ∇∇:
        # https://github.com/gridap/Gridap.jl/blob/7f5f15e53cd390b6ac82d51ee46fa5bc0792a7e3/src/Fields/AutoDiff.jl#L9
        # Definition of  Δ:
        # https://github.com/gridap/Gridap.jl/blob/3916f3c87d86dfc6deb2b24be0305614267c5785/src/Fields/DiffOperators.jl#L67

        # PROBLEM 2
        # Statement: Following up from the previous problem. It should be corrected to
        # mean(Δ(u)) -> mean(n_Λ*∇∇(u)*n_Λ)

        # PROBLEM 3
        # Statement: Why is the terms negative??

        a(u,v) = ∫( Δ(u)*Δ(v) + α* u⋅v )dΩ +
                 ∫( - mean(Δ(u))*jump(∇(v)⋅n_Λ) - jump(∇(u)⋅n_Λ)*mean(Δ(v))
                   + γ/h*jump(∇(u)⋅n_Λ)*jump(∇(v)⋅n_Λ) )dΛ

        # PROBLEM 4
        # Why the directional derivative of test function v? Does not makes sense given the identity:
        # --> (Δ^2 u, v)_Ω  = (D^2 u , D^2 v)_Ω + (∂_n Δ u, v)_∂Ω - (∂_nn u, ∂_n v)_∂Ω  - (∂_nt u, ∂_t v)_∂Ω
        #                   = (D^2 u , D^2 v)_Ω + (g, v)_∂Ω

        l(v) = ∫( v*f )dΩ + ∫( g*(∇(v)⋅n_Γ) )dΓ
        op = AffineFEOperator(a,l,U,V)
        uh = solve(op)

        u_inter = interpolate_everywhere(u, V)

        e = u_inter - uh
        l2(u) = sqrt(sum( ∫( u⊙u )*dΩ ))
        h_energy(u) = sqrt(sum( ∫( ∇(u)⊙∇(u) )*dΩ
                         + ( γ/h ) * ∫(jump(∇(e)⋅n_Λ) ⊙ jump(∇(e)⋅n_Λ))dΛ
                         + ( h/γ ) * ∫(mean(Δ(e)) ⊙ mean(Δ(e)))dΛ
                        ))

        el2 = l2(e)
        eh = h_energy(e)
        eh = el2
        res = Results( model=model, Ω=Ω, Γ=Γ, Λ=Λ,
                      h=h, γ=γ, order=order, degree=degree,
                      u_inter=u_inter, uh=uh, e=e, el2=el2, eh=eh)

        return res
    end


    function man_sol(;L=1,m=1,r=1)
        u(x) = cos(m*( 2π/L )*x[1])*cos(r*( 2π/L )*x[2])
    end

end # module






