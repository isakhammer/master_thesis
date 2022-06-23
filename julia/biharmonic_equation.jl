

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
        eh1
    end


    function generate_square_spaces(;n, L, γ=1, order=2, simplex=true, u=nothing)
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

        return BiharmonicEquation.GridapSpaces(model, h, γ, order, degree,
                                               V,U,
                                               Ω, Γ, Λ,
                                               dΩ, dΓ, dΛ,
                                               n_Γ, n_Λ)
    end

    function generate_sol(;u, uh, ss::GridapSpaces)
        e = u - uh
        l2(u) = sqrt(sum( ∫( u⊙u )*ss.dΩ ))
        h1(u) = sqrt(sum( ∫( u⊙u + ∇(u)⊙∇(u) )*ss.dΩ ))
        el2 = l2(e)
        eh1 = h1(e)
        Solution(u,uh,e,el2,eh1)
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

        fig = plot(ss.Ω)
        wireframe!(ss.Ω, color=:black, linewidth=2)
        # scatter!(ss.Ω, marker=:star8, markersize=20, color=:blue)
        save(dirname*"/grid.png", fig)

        fig = plot(ss.Λ)
        wireframe!(ss.Λ, color=:black, linewidth=2)
        # scatter!(ss.Ω, marker=:star8, markersize=20, color=:blue)
        save(dirname*"/lambda.png", fig)

        fig, _ , plt = plot(ss.Ω, sol.u)
        Colorbar(fig[1,2], plt)
        save(dirname*"/man_sol.png", fig)

    end


    function run_test_method(;ss::GridapSpaces, u::Function)
        # Analytical manufactured solution
        α = 1

        f(x) = Δ(Δ(u))(x)+ α*u(x)
        g(x) = Δ(u)(x)

        # @test f(VectorValue(0.5,0.5)) == ( 4+α )*u(VectorValue(0.5,0.5)) # redo rhs
        # @test g(VectorValue(0.5,0.5)) == -2*u(VectorValue(0.5,0.5))      # redo rhs

        # Weak form
        γ = ss.γ

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
                   + γ/ss.h*jump(∇(u)⋅ss.n_Λ)*jump(∇(v)⋅ss.n_Λ) )ss.dΛ

        # PROBLEM 3
        # Why the directional derivative of test function v? Does not makes sense given the identity
        # (Δ^2 u, v)_Ω = (D^2 u , D^2 v)_Ω + (∂_n Δ u, v)_∂Ω - (∂_nn u, ∂_n v)_∂Ω  - (∂_nt u, ∂_t v)_∂Ω
        #              = (D^2 u , D^2 v)_Ω + (g, v)_∂Ω

        l(v) = ∫( v*f )ss.dΩ + ∫( g*(∇(v)⋅ss.n_Γ) )ss.dΓ
        op = AffineFEOperator(a,l,ss.U,ss.V)

        uh = solve(op)
        sol = generate_sol(u=u,uh=uh,ss=ss)
        return sol
    end

    function run_DG_method(;ss::GridapSpaces, u::Function)
        f(x) = Δ(Δ(u))(x)+ α*u(x)
        g(x) = Δ(u)(x)

        α = 1

        mean_nn(u) = 0.5*( ss.n_Λ.plus⋅ ∇∇(u).plus⋅ ss.n_Λ.plus + ss.n_Λ.minus ⋅ ∇∇(u).minus ⋅ ss.n_Λ.minus )
        jump_n(u) = ∇(u).plus⋅ ss.n_Λ.plus - ∇(u).minus ⋅ ss.n_Λ.minus
⋅
        # Inner facets
        a(u,v) = ∫( ∇∇(v)⊙∇∇(u) + α⋅(v⊙u) )ss.dΩ + ∫(  mean_nn(v)⊙jump_n(u) + mean_nn(u)⊙jump_n(v) + (ss.γ/ss.h)⋅ jump_n(v)⊙jump_n(u))ss.dΛ

        l(v) = ∫( v ⋅ f )ss.dΩ + ∫(- (g⋅v))ss.dΓ

        op = AffineFEOperator(a, l, ss.U, ss.V)
        uh = solve(op)

        sol = generate_sol(u=u,uh=uh,ss=ss)
        return sol
    end

    function man_sol(;L=1,m=1,r=1)
        u(x) = cos(m*( 2π/L )*x[1])*cos(r*( 2π/L )*x[2])
    end



    function run_CP_method(;ss::GridapSpaces, u::Function, method="test")
        method=="test" && return run_test_method(ss=ss, u=u)
        method=="DG" && return run_DG_method(ss=ss, u=u)
        throw(DomainError(method, "Does not have a method with this name"))
    end

end # module






