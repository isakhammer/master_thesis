using Gridap
using Plots
using Test

function run_poission(; n=10, generate_vtk::Bool=true, dirname::String, test::Bool = false)

    # source: https://gridap.github.io/Tutorials/v0.5/pages/t005_dg_discretization/

    # Problem statement:
    # -Δu = f
    # u = g
    # Ω = (0,1)^2

    u(x) = 3*x[1] + x[2]^2
    f(x) = - Δ(u)(x)
    g(x) = u(x)

    # mesh generation
    L = 1.0
    order = 1

    domain2D = (0.0, L, 0.0, L)
    partition2D = (n,n)
    model = CartesianDiscreteModel(domain2D,partition2D)

    V = TestFESpace(model, ReferenceFE(lagrangian,Float64,order), conformity=:L2)
    U = TrialFESpace(V)

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
    a_Ω(u,v) = ∫( ∇(v)⊙∇(u) )dΩ
    l_Ω(v) = ∫( v*f )dΩ

    h = L / n
    γ = order*(order+1)

    a_Γ(u,v) = ∫( - v*(∇(u)⋅n_Γ) - (∇(v)⋅n_Γ) * u + (γ/h)*v*u )dΓ
    l_Γ(v)   = ∫(                - (∇(v)⋅n_Γ)*g + (γ/h)*v*g )dΓ


    a_Λ(u,v) = ∫( - jump(v*n_Λ)⊙mean(∇(u))
              - mean(∇(v))⊙jump(u*n_Λ)
              + (γ/h)*jump(v*n_Λ)⊙jump(u*n_Λ) )dΛ

    a(u,v) = a_Ω(u,v) + a_Γ(u,v) + a_Λ(u,v)
    l(v) = l_Ω(v) + l_Γ(v)

    op = AffineFEOperator(a, l, U, V)
    uh = solve(op)

    e = u - uh
    el2 = sqrt(sum( ∫(e*e)dΩ ))
    eh1 = sqrt(sum( ∫( e*e + ∇(e)⋅∇(e) )*dΩ ))

    if test==true
        @test el2 < 10^-3
    end

    if !generate_vtk
        return el2, eh1
    end


    # Generate plots
    writevtk(model, dirname*"/model")
    writevtk(Λ,dirname*"/skeleton")
    writevtk(Λ,dirname*"/jumps",cellfields=["jump_u"=>jump(uh)])
    writevtk(Ω,dirname*"/omega",cellfields=["uh"=>uh])
    writevtk(Ω,dirname*"/error",cellfields=["e"=>e])
    writevtk(Ω,dirname*"/manufatured",cellfields=["u"=>u])
    return el2, eh1

end

function conv_test(;dirname)
    ns = [8,16,32,64,128]

    el2s = Float64[]
    eh1s = Float64[]
    hs = Float64[]

    println("Run convergence tests")
    for n in ns

        el2, eh1 = run_poission(n=n,dirname=dirname)
        println("Simulation with n:", n, ", Errors:  L2: ", el2, " H1:", eh1)

        h = ( 1/n )

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
    dirname = "poission_results"

    if (isdir(dirname))
        rm(dirname, recursive=true)
    end
    mkdir(dirname)

    run_poission(n=20, generate_vtk=true, dirname=dirname, test=true)
    conv_test(dirname=dirname)
end

main()

