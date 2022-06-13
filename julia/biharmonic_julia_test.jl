
using Test
using Gridap
using Plots
gr()

function run_biharmonic_julia_test(; n=10, generate_vtk=false, dirname="biharmonic_results", test=false)
    # Analytical manufactured solution
    α = 1

    # u(x) = x[1]*(x[1]-1)*x[2]*(x[2]-1)
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
    nΓ = get_normal_vector(Γ)
    nΛ = get_normal_vector(Λ)

    # Weak form
    h = (domain[2]-domain[1]) / partition[1]
    γ = 1
    a(u,v) = ∫( Δ(u)*Δ(v) + α* u⋅v )dΩ +
             ∫( - mean(Δ(u))*jump(∇(v)⋅nΛ) - jump(∇(u)⋅nΛ)*mean(Δ(v)) + γ/h*jump(∇(u)⋅nΛ)*jump(∇(v)⋅nΛ) )dΛ
    l(v) = ∫( v*f )dΩ + ∫( g*(∇(v)⋅nΓ) )dΓ
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
