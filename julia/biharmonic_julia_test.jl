

using Test
using Gridap
using Plots
using LaTeXStrings
# plotlyjs()

module BiharmonicEquation
    using Gridap

    struct GridapSpaces
        model
        h::Float64
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

    function generate_square_space(;n, L=2π, order=2, simplex=true)
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
        # U = TrialFESpace(V,u)
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

        return BiharmonicEquation.GridapSpaces(model, h, order, degree,
                                               V,U,
                                               Ω, Γ, Λ,
                                               dΩ, dΓ, dΛ,
                                               n_Γ, n_Λ)
    end

    struct Solution
        u
        uh
        e
        el2
        eh1
    end

    function generate_sol(u,uh, ss::GridapSpaces)
        e = u - uh
        l2(u) = sqrt(sum( ∫( u⊙u )*ss.dΩ ))
        h1(u) = sqrt(sum( ∫( u⊙u + ∇(u)⊙∇(u) )*ss.dΩ ))
        el2 = l2(e)
        eh1 = h1(e)
        Solution(u,uh,e,el2,eh1)
    end

end




function run_biharmonic_julia_test(; n=10, order::Int, generate_vtk=false,
        dirname="biharmonic_julia_test_results", test=false, simplex=false)
    # Analytical manufactured solution
    α = 1
    u(x) = cos(x[1])*cos(x[2])
    f(x) = Δ(Δ(u))(x)+ α*u(x)
    g(x) = Δ(u)(x)


    # Domain
    @test f(VectorValue(0.5,0.5)) == ( 4+α )*u(VectorValue(0.5,0.5))
    @test g(VectorValue(0.5,0.5)) == -2*u(VectorValue(0.5,0.5))


    ss = BiharmonicEquation.generate_square_space(n=n, order=order)

    # Weak form
    # h = (domain[2]-domain[1]) / partition[1]
    γ = 1

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

    sol = BiharmonicEquation.generate_sol(u,uh,ss)

    # Error
    # tol = 1.0e-10
    # @test el2 < tol
    # @test eh1 < tol

    if test==true
        @test sol.el2 < 10^-1
    end

    if !generate_vtk
        return sol.el2, sol.eh1
    end

    writevtk(ss.model, dirname*"/model")
    writevtk(ss.Λ,dirname*"/skeleton")
    writevtk(ss.Γ,dirname*"/boundary")
    writevtk(ss.Λ,dirname*"/jumps",cellfields=["jump_u"=>jump(sol.uh)])
    writevtk(ss.Ω,dirname*"/omega",cellfields=["uh"=>sol.uh])
    writevtk(ss.Ω,dirname*"/error",cellfields=["e"=>sol.e])
    writevtk(ss.Ω,dirname*"/manufatured",cellfields=["u"=>sol.u])

    return sol.el2, sol.eh1

end

function conv_test(; dirname, order)
    ns = [8,16,25,32,64]
    # ns = collect(10:10:100)

    el2s = Float64[]
    eh1s = Float64[]
    hs = Float64[]

    println()
    println("Run convergence tests: order = "*string(order))

    for n in ns
        el2, eh1 = run_biharmonic_julia_test(n=n, order=order)
        println("Simulation with n:", n, ", Errors:  L2: ", el2, " H1:", eh1)
        h = ( 1/n )*2*π
        push!(el2s,el2)
        push!(eh1s,eh1)
        push!(hs,h)

    end

    function slope(hs,errors)
      x = log10.(hs)
      y = log10.(errors)
      linreg = hcat(fill!(similar(x), 1), x) \ y
      linreg[2]
    end

    p_L2 = slope(hs,el2s)
    p_H1 = slope(hs,eh1s)

    println("Slope of L2 is ", p_L2)
    println("Slope of H1 is ", p_H1)

    p = Plots.plot(hs,[el2s eh1s ],
        xaxis=:log, yaxis=:log,
        label=[L"Error norm in $L_2(\Omega)$ where $p_1 = $"*string(round(p_L2,digits=2)) L"Error norm in  $ H^1(\Omega)$ where $p_2 =$"*string(round(p_H1,digits=2)) ],
        shape=:auto,
        legend=:topleft,
        legendfontsize=10,
        xlabel=L"$h$",ylabel="error norm" , show = true)

    Plots.savefig(p, dirname*"/convergence_d_"*string(order)*".png")
end



function main()

    # Generate plots
    function makedir(dirname)
        if (isdir(dirname))
            rm(dirname, recursive=true)
        end
        mkdir(dirname)
    end

    folder = "biharmonic_julia_test_results"
    makedir(folder)

    exampledir = folder*"/example"
    makedir(exampledir)

    println("Generating examples")
    ns = [100]
    for n in ns
        ndir = exampledir*"/n_"*string(n)
        makedir(ndir)
        run_biharmonic_julia_test(n=n, order=2, generate_vtk=true, dirname=ndir, test=true,simplex=true)
    end

    println("Generating convergence tests")
    plotdir = folder*"/plots"
    makedir(plotdir)
    # orders = [1,2,3]
    orders = [1,2,3,4]
    for order in orders
        conv_test(dirname=plotdir, order=order)
    end
end



main()
