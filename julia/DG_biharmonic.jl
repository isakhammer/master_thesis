using Gridap
import Gridap: ∇

function main()


    # mesh generation
    L = 2*π
    n = 13
    α = 1
    h = L / n
    γ = 2
    order = 2
    domain2D = (0, L, 0, L)
    partition2D = (n,n)
    model = CartesianDiscreteModel(domain2D,partition2D)
    writevtk(model,"plots/biharmonic_model")


    # Spaces
    V = TestFESpace(model, ReferenceFE(lagrangian,Float64,order), conformity=:H1)
    U = TrialFESpace(V)
    Ω = Triangulation(model)
    Γ = BoundaryTriangulation(model)
    Λ = SkeletonTriangulation(model)

    writevtk(Λ,"plots/biharmonic_skeleton")

    degree = 1*order

    dΩ = Measure(Ω,degree)
    dΓ = Measure(Γ,degree)
    dΛ = Measure(Λ,degree)

    n_Γ = get_normal_vector(Γ)
    n_Λ = get_normal_vector(Λ)

    # manufactured solution
    u(x) = cos(x[1])*cos(x[2])

    # f(x) = (α + 4)* cos(x[1])*cos(x[2])
    f(x) = Δ(u)(x) + α*u(x) # Algorithmic Diff.
    g(x) = 1  # we see that u_h -> 0 when g -> 0

    # Inner triangulation
    a_Ω(u,v) = ∫( ∇∇(v)⊙∇∇(u) + α⋅(v⊙u) )dΩ
    l_Ω(v) = ∫( v ⋅ f )dΩ

    # Outer facets
    l_Γ(v) = ∫(- (g⋅v))dΓ

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
    l(v) = l_Ω(v) + l_Γ(v)

    op = AffineFEOperator(a, l, U, V)
    uh = solve(op)
    writevtk(Λ,"plots/biharmonic_jumps",cellfields=["jump_u"=>jump(uh)])
    writevtk(Ω,"plots/biharmonic_omega",cellfields=["uh"=>uh])
    writevtk(Ω,"plots/biharmonic_manufatured",cellfields=["u"=>u])
end

main()
