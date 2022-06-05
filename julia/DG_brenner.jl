using Gridap
using Plots
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
    U = TrialFESpace(V, 0.0)
    Ω = Triangulation(model)
    Γ = BoundaryTriangulation(model)
    Λ = SkeletonTriangulation(model)

    degree = 1*order
    dΩ = Measure(Ω,degree)
    dΓ = Measure(Γ,degree)
    dΛ = Measure(Λ,degree)

    n_Λ = get_normal_vector(Λ)

    f=1

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

    if !generate_vtk
        return 0
    end

    if (isdir(dirname))
        rm(dirname, recursive=true)
    end

    mkdir(dirname)
    writevtk(model, dirname*"/biharmonic_model")
    writevtk(Λ,dirname*"/brenner_skeleton")
    writevtk(Λ,dirname*"/brenner_jumps",cellfields=["jump_u"=>jump(uh)])
    writevtk(Ω,dirname*"/brenner_omega",cellfields=["uh"=>uh])

    return
end

run_brenner(n=40, generate_vtk=true, dirname="brenner_results", test=false)

