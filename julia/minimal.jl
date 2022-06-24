using Gridap
using Plots
gr()
using Test
import Gridap: ∇

function run_CP(; n=10, generate_vtk::Bool=false, dirname::String, test::Bool=false)

    # mesh generation
    L = 2π
    h = L / n
    γ = 0.5
    u(x) = cos(x[1])*cos(x[2])

    order = 2
    domain2D = (0, L, 0, L)
    partition2D = (n,n)
    model = CartesianDiscreteModel(domain2D,partition2D)

    # Spaces
    V = TestFESpace(model, ReferenceFE(lagrangian,Float64,order), conformity=:L2)
    U = TrialFESpace(V)
    Ω = Triangulation(model)
    Γ = BoundaryTriangulation(model)
    Λ = SkeletonTriangulation(model)

    degree = 2*order
    dΩ = Measure(Ω,degree)
    dΓ = Measure(Γ,degree)
    dΛ = Measure(Λ,degree)

    n_Λ = get_normal_vector(Λ)

    # manufactured solution
    f(x) = Δ(Δ(u))(x)+ α*u(x)
    g(x) = Δ(u)(x)
    α = 1

    mean_nn(u,n) = mean(n ⋅ (∇∇(u) ⋅ n))
    mean_nn(u,n) = 0.5*( n.plus⋅ ∇∇(u).plus⋅ n.plus + n.minus ⋅ ∇∇(u).minus ⋅ n.minus )
    ⋅
    # Inner facets
    a(u,v) = ∫( ∇∇(v)⊙∇∇(u) + α⋅(v⊙u) )dΩ + ∫(-mean_nn(v,n_Λ)⊙jump(∇(u)⋅n_Λ) - mean_nn(u,n_Λ)⊙jump(∇(v)⋅n_Λ))dΛ + ∫((γ/h)⋅jump(∇(u)⋅n_Λ)⊙jump(∇(v)⋅n_Λ))dΛ
        # + ∫(-mean_nn(v,n_Γ )⊙jump(∇(u)⋅n_Γ) - mean_nn(u,n_Γ)⊙jump(∇(v)⋅n_Γ))dΓ + ∫((γ/h)⋅jump(∇(u)⋅n_Γ)⊙jump(∇(v)⋅n_Γ))dΓ

    l(v) = ∫( v ⋅ f )dΩ + ∫(- (g⋅v))dΓ


    op = AffineFEOperator(a, l, U, V)
    uh = solve(op)

    e = u - uh
    el2 = sqrt(sum( ∫(e*e)dΩ ))
    eh1 = sqrt(sum( ∫( e*e + ∇(e)⋅∇(e) )*dΩ ))

    if !generate_vtk
        return el2, eh1
    end

    writevtk(model, dirname*"/model")
    writevtk(Λ,dirname*"/skeleton")
    writevtk(Λ,dirname*"/jumps",cellfields=["jump_u"=>jump(uh)])
    writevtk(Ω,dirname*"/omega",cellfields=["uh"=>uh])
    writevtk(Ω,dirname*"/error",cellfields=["e"=>e])
    writevtk(Ω,dirname*"/manufatured",cellfields=["u"=>u])

    if test==true
        @test el2 < 10^-5
    end
    return el2, eh1
end



function conv_test(;dirname)
    ns = [8,16,32,64,128]

    el2s = Float64[]
    eh1s = Float64[]
    hs = Float64[]

    println("Run convergence tests")
    for n in ns

        el2, eh1 = run_CP(n=n, dirname=dirname)
        println("Simulation with n:", n, ", Errors:  L2: ", el2, " H1:", eh1)

        h = ( 1/n )*2*π
        push!(el2s,el2)
        push!(eh1s,eh1)
        push!(hs,h)
    end

    p = Plots.plot(hs,[el2s eh1s],
                      # xaxis=:log, yaxis=:log,
                      label=["L2" "H1"],
                      shape=:auto,
                      xlabel="h",ylabel="error norm")

    Plots.savefig(p, dirname*"/convergence.png")
end

function main()
    dirname = "minimal_example"

    if (isdir(dirname))
        rm(dirname, recursive=true)
    end
    mkdir(dirname)

    # run_CP(n=10, generate_vtk=true, dirname=dirname, test=false)
    conv_test(dirname=dirname)
end

main()
