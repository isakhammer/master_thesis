
using Test
using Gridap
using Plots
# plotlyjs()


function run_biharmonic_julia_test(; n=10, generate_vtk=false, dirname="biharmonic_julia_test_results", test=false)
    # Analytical manufactured solution
    α = 1

    u(x) = cos(x[1])*cos(x[2])

    f(x) = Δ(Δ(u))(x)+ α*u(x)
    g(x) = Δ(u)(x)


    # Domain
    L = 2*π

    @test f(VectorValue(0.5,0.5)) == ( 4+α )*u(VectorValue(0.5,0.5))
    @test g(VectorValue(0.5,0.5)) == -2*u(VectorValue(0.5,0.5))

    domain = (0,L,0,L)
    partition = (n,n)
    model = CartesianDiscreteModel(domain,partition)

    # FE space
    order = 2
    V = TestFESpace(model,ReferenceFE(lagrangian,Float64,order),dirichlet_tags="boundary")
    U = TrialFESpace(V,u)

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

    # Weak form
    h = (domain[2]-domain[1]) / partition[1]
    h = L / n
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

    a(u,v) = ∫( Δ(u)*Δ(v) + α* u⋅v )dΩ +
             ∫( - mean(Δ(u))*jump(∇(v)⋅n_Λ) - jump(∇(u)⋅n_Λ)*mean(Δ(v))
               + γ/h*jump(∇(u)⋅n_Λ)*jump(∇(v)⋅n_Λ) )dΛ

    # PROBLEM 3
    # Why the directional derivative of test function v? Does not makes sense given the identity
    # (Δ^2 u, v)_Ω = (D^2 u , D^2 v)_Ω + (∂_n Δ u, v)_∂Ω - (∂_nn u, ∂_n v)_∂Ω  - (∂_nt u, ∂_t v)_∂Ω
    #              = (D^2 u , D^2 v)_Ω + (g, v)_∂Ω

    l(v) = ∫( v*f )dΩ + ∫( g*(∇(v)⋅n_Γ) )dΓ
    op = AffineFEOperator(a,l,U,V)

    uh = solve(op)

    # Error
    e = u - uh
    l2(u) = sqrt(sum( ∫( u⊙u )*dΩ ))
    h1(u) = sqrt(sum( ∫( u⊙u + ∇(u)⊙∇(u) )*dΩ ))
    el2 = l2(e)
    eh1 = h1(e)
    # tol = 1.0e-10
    # @test el2 < tol
    # @test eh1 < tol

    if test==true
        @test el2 < 10^-3
    end

    if !generate_vtk
        return el2, eh1
    end

    writevtk(model, dirname*"/model")
    writevtk(Λ,dirname*"/skeleton")
    writevtk(Λ,dirname*"/jumps",cellfields=["jump_u"=>jump(uh)])
    writevtk(Ω,dirname*"/omega",cellfields=["uh"=>uh])
    writevtk(Ω,dirname*"/error",cellfields=["e"=>e])
    writevtk(Ω,dirname*"/manufatured",cellfields=["u"=>u])

    return el2, eh1

end

function conv_test(; dirname)
    ns = [8,16,32,64,128]

    el2s = Float64[]
    eh1s = Float64[]
    hs = Float64[]

    println("Run convergence tests")
    for n in ns

        el2, eh1 = run_biharmonic_julia_test(n=n)
        println("Simulation with n:", n, ", Errors:  L2: ", el2, " H1:", eh1)

        h = ( 1/n )*2*π

        push!(el2s,el2)
        push!(eh1s,eh1)
        push!(hs,h)

    end

    p = Plots.plot(hs,[el2s eh1s],
        xaxis=:log, yaxis=:log,
        label=["L2" "H1"],
        shape=:auto,
        xlabel="h",ylabel="error norm")

    function slope(hs,errors)
      x = log10.(hs)
      y = log10.(errors)
      linreg = hcat(fill!(similar(x), 1), x) \ y
      linreg[2]
    end

    println("Slope of L2 is ", slope(hs,el2s))
    println("Slope of H1 is ", slope(hs,eh1s))
    Plots.savefig(p, dirname*"/convergence.png")
end


function main()
    dirname = "biharmonic_julia_test_results"

    # Generate plots
    if (isdir(dirname))
        rm(dirname, recursive=true)
    end
    mkdir(dirname)

    run_biharmonic_julia_test(n=90, generate_vtk=true, dirname=dirname, test=true)
    conv_test(dirname=dirname)
end

main()
