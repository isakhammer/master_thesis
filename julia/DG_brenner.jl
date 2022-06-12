using Gridap

using Plots
gr()
using Test
import Gridap: ∇

function run_brenner(; n=10, generate_vtk=false, dirname="biharmonic_results", test=false)

    # mesh generation
    L = 1
    α = 1
    h = L / n
    γ = 1

    order = 2
    domain2D = (0, L, 0, L)
    partition2D = (n,n)
    model = CartesianDiscreteModel(domain2D,partition2D)

    # Spaces
    V = TestFESpace(model, ReferenceFE(lagrangian,Float64,order), conformity=:H1, dirichlet_tags= "boundary")
    dirichlet_condition = 0
    U = TrialFESpace(V, dirichlet_condition)
    Ω = Triangulation(model)
    Γ = BoundaryTriangulation(model)
    Λ = SkeletonTriangulation(model)

    degree = 1*order
    dΩ = Measure(Ω,degree)
    dΓ = Measure(Γ,degree)
    dΛ = Measure(Λ,degree)

    n_Λ = get_normal_vector(Λ)

    # manufactured solution
    u(x) = 400*( x[1]*(x[1]-1) )^2*( x[2]*(x[2]-1) )^2
    f(x) = Δ(Δ(u))(x)
    # f(x) = 2*x[1]^2 *(x[1] - 1)^2*(x[2]^2 + 4*x[2]*(x[2] - 1) + (x[2] - 1)^2) + 2*x[2]^2*(x[2] - 1)^2*(x[1]^2 + 4*x[1]*(x[1] - 1) + (x[1] - 1)^2)



    # Inner triangulation
    a_Ω(u,v) = ∫( ∇∇(v)⊙∇∇(u) + α* u ⊙ v  )dΩ
    l_Ω(v) = ∫( v ⋅ f )dΩ

    mean_nn(u) = 0.5*( n_Λ.plus⋅ ∇∇(u).plus⋅ n_Λ.plus + n_Λ.minus ⋅ ∇∇(u).minus ⋅ n_Λ.minus )
    jump_n(u) = ∇(u).plus⋅ n_Λ.plus - ∇(u).minus ⋅ n_Λ.minus
⋅
    # Inner facets
    a_Λ(u,v) = ∫(  mean_nn(v)⊙jump_n(u)
                 + mean_nn(u)⊙jump_n(v)
                 + (γ/h)⋅ jump_n(v)⊙jump_n(u)
                )dΛ

    # Summation
    a(u,v) = a_Ω(u,v) + a_Λ(u,v)
    l(v) = l_Ω(v)

    op = AffineFEOperator(a, l, U, V)
    uh = solve(op)

    e = u - uh
    el2 = sqrt(sum( ∫(e*e)dΩ ))
    eh1 = sqrt(sum( ∫( e*e + ∇(e)⋅∇(e) )*dΩ ))



    if !generate_vtk
        return el2, eh1
    end

    if (isdir(dirname))
        rm(dirname, recursive=true)
    end

    mkdir(dirname)
    writevtk(model, dirname*"/brenner_model")
    writevtk(Λ,dirname*"/brenner_skeleton")
    writevtk(Λ,dirname*"/brenner_jumps",cellfields=["jump_u"=>jump(uh)])
    writevtk(Ω,dirname*"/brenner_omega",cellfields=["uh"=>uh])
    writevtk(Ω,dirname*"/brenner_error",cellfields=["e"=>e])
    writevtk(Ω,dirname*"/brenner_manufatured",cellfields=["u"=>u])

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

        el2, eh1 = run_brenner(n=n)
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
    dirname = "brenner_results"

    if (isdir(dirname))
        rm(dirname, recursive=true)
    end
    mkdir(dirname)

    run_biharmonic(n=10, generate_vtk=true, dirname=dirname, test=false)
    conv_test(dirname=dirname)
end

main()
