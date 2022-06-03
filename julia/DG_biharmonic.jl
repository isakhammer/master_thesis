using Gridap
import Gridap: ∇

function main()
    f(x) = 1
    g(x) = 1

    # mesh generation
    L = 1.0
    n = 5
    α = 3
    order = 2

    domain2D = (0.0, L, 0.0, L)
    partition2D = (n,n)
    model = CartesianDiscreteModel(domain2D,partition2D)
    writevtk(model,"plots/biharmonic_model")

    V = TestFESpace(model, ReferenceFE(lagrangian,Float64,order), conformity=:H1)
    U = TrialFESpace(V)
    Ω = Triangulation(model)
    Γ = BoundaryTriangulation(model)
    Λ = SkeletonTriangulation(model)

    writevtk(Λ,"plots/biharmonic_skeleton")

    degree = 2*order

    dΩ = Measure(Ω,degree)
    dΓ = Measure(Γ,degree)
    dΛ = Measure(Λ,degree)

    n_Γ = get_normal_vector(Γ)
    n_Λ = get_normal_vector(Λ)

    # Inner triangulation
    a_Ω(u,v) = ∫( ∇∇(v)⊙∇∇(u) + α⋅(v⊙u) )dΩ
    l_Ω(v) = ∫( v ⋅ f )dΩ

    # Outer facets
    l_Γ(v) = ∫(- (g⋅v))dΓ

    h = L / n
    γ = 2

    mean_∇∇(u) = 0.5*( n_Λ.plus⋅ ∇∇(u).plus⋅ n_Λ.plus + n_Λ.minus ⋅ ∇∇(u).minus ⋅ n_Λ.minus )
    mean_∇(u) = 0.5*(  ∇(u).plus⋅ n_Λ.plus + ∇(u).minus ⋅ n_Λ.minus )
    jump_∇(u) = ∇(u).plus⋅ n_Λ.plus - ∇(u).minus ⋅ n_Λ.minus
⋅
    # Inner facets
    a_Λ(u,v) = ∫(
                 + mean_∇∇(v)⊙jump_∇(u)
                 + mean_∇∇(u)⊙jump_∇(v)
              + (γ/h)⋅ jump_∇(v)⊙jump_∇(u)
                )dΛ

    # Summation
    a(u,v) = a_Ω(u,v) + a_Λ(u,v)
    l(v) = l_Ω(v) + l_Γ(v)

    op = AffineFEOperator(a, l, U, V)
    uh = solve(op)
    writevtk(Λ,"plots/biharmonic_jumps",cellfields=["jump_u"=>jump(uh)])
    writevtk(Ω,"plots/biharmonic_omega",cellfields=["uh"=>uh])
end

main()
