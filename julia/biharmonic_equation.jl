

module BiharmonicEquation
    using Gridap
    using GridapMakie
    using GLMakie
    using Test

    # TODO: find a way to make structs in a one-liner
    struct GridapSpaces
        model
        h::Real
        γ::Real
        order::Int
        degree::Int

        V
        U

        Ω
        Γ
        Λ

        dΩ
        dΓ
        dΛ

        n_Γ
        n_Λ
    end

    struct Solution
        u
        uh
        e
        el2
        eh
    end


    function generate_vtk(;ss::GridapSpaces, sol::Solution, dirname::String)
        println("Generating vtk's in ", dirname)
        if (isdir(dirname))
            rm(dirname, recursive=true)
        end
        mkdir(dirname)

        writevtk(ss.model, dirname*"/model")
        writevtk(ss.Λ,dirname*"/skeleton")
        writevtk(ss.Γ,dirname*"/boundary")
        writevtk(ss.Λ,dirname*"/jumps",cellfields=["jump_u"=>jump(sol.uh)])
        writevtk(ss.Ω,dirname*"/omega",cellfields=["uh"=>sol.uh])
        writevtk(ss.Ω,dirname*"/error",cellfields=["e"=>sol.e])
        writevtk(ss.Ω,dirname*"/manufatured",cellfields=["u"=>sol.u])


        fig = plot(ss.Λ)
        wireframe!(ss.Λ, color=:black, linewidth=2)
        wireframe!(ss.Γ, color=:black, linewidth=2)
        save(dirname*"/grid.png", fig)

        fig, _ , plt = plot(ss.Ω, sol.u)
        Colorbar(fig[1,2], plt)
        save(dirname*"/man_sol.png", fig)

    end


    function run_CP_method(;n, L, γ, order,  u::Function = man_sol(), simplex=true,)
        h = L / n

        domain = (0,L,0,L)
        partition = (n,n)

        if simplex
            model = CartesianDiscreteModel(domain,partition) |> simplexify
        else
            model = CartesianDiscreteModel(domain,partition)
        end

        # FE space
        V = TestFESpace(model,ReferenceFE(lagrangian,Float64,order))

        if u !== nothing
            U = TrialFESpace(V,u)
        else
            U = TrialFESpace(V)
        end

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

        ss = BiharmonicEquation.GridapSpaces(model, h, γ, order, degree,
                                                    V,U,
                                                    Ω, Γ, Λ,
                                                    dΩ, dΓ, dΛ,
                                                    n_Γ, n_Λ)

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

        a(u,v) = ∫( Δ(u)*Δ(v) + α* u⋅v )ss.dΩ +
                 ∫( - mean(Δ(u))*jump(∇(v)⋅ss.n_Λ) - jump(∇(u)⋅ss.n_Λ)*mean(Δ(v))
                   + ss.γ/ss.h*jump(∇(u)⋅ss.n_Λ)*jump(∇(v)⋅ss.n_Λ) )ss.dΛ

        # PROBLEM 4
        # Why the directional derivative of test function v? Does not makes sense given the identity:
        # --> (Δ^2 u, v)_Ω  = (D^2 u , D^2 v)_Ω + (∂_n Δ u, v)_∂Ω - (∂_nn u, ∂_n v)_∂Ω  - (∂_nt u, ∂_t v)_∂Ω
        #                   = (D^2 u , D^2 v)_Ω + (g, v)_∂Ω

        l(v) = ∫( v*f )ss.dΩ + ∫( g*(∇(v)⋅ss.n_Γ) )ss.dΓ
        op = AffineFEOperator(a,l,ss.U,ss.V)
        uh = solve(op)

        e = u - uh
        l2(u) = sqrt(sum( ∫( u⊙u )*ss.dΩ ))
        h_energy(u) = sqrt(sum( ∫( ∇(u)⊙∇(u) )*ss.dΩ
                         + ( ss.γ/ss.h ) * ∫(jump(∇(e)⋅ss.n_Λ) ⊙ jump(∇(e)⋅ss.n_Λ))ss.dΛ
                         + ( ss.h/ss.γ ) * ∫(mean(Δ(e)) ⊙ mean(Δ(e)))ss.dΛ
                        ))

        el2 = l2(e)
        eh = h_energy(e)
        sol = Solution(u,uh,e,el2,eh)
        return sol, ss
    end


    function man_sol(;L=1,m=1,r=1)
        u(x) = cos(m*( 2π/L )*x[1])*cos(r*( 2π/L )*x[2])
    end

end # module






