using Gridap
import Gridap: ∇



function main()
    f(x) = 1
    g(x) = 0

    # ∇u(x) = VectorValue(3   + x[2]*x[3],
    #                 2*x[2]   + x[1]*x[3],
    #                 6*x[3]^2 + x[1]*x[2])

    # println(∇(u) === ∇u)

    # mesh generation
    L = 1.0
    n = 4
    α = 1
    order = 2

    domain2D = (0.0, L, 0.0, L)
    partition2D = (n,n)
    model = CartesianDiscreteModel(domain2D,partition2D)
    writevtk(model,"plots/model_mesh")

    V = TestFESpace(model, ReferenceFE(lagrangian,Float64,order), conformity=:L2)
    U = TrialFESpace(V)
    Ω = Triangulation(model)
    Γ = BoundaryTriangulation(model)

    Λ = SkeletonTriangulation(model)

    writevtk(Λ,"plots/strian_cube")

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
    γ = order*(order+1)

    # Inner facets
    a_Λ(u,v) = ∫(
                 + mean(n_Λ⋅ ∇∇(v)⋅ n_Λ)⊙jump(∇(u)⋅n_Λ)
                 + mean(n_Λ⋅ ∇∇(u)⋅ n_Λ)⊙jump(∇(v)⋅n_Λ)
              + (γ/h)⋅ mean(∇(v) ⋅ n_Λ)⊙jump(∇(u)⋅ n_Λ)
                )dΛ

    # Summation
    a(u,v) = a_Ω(u,v) + a_Λ(u,v)
    l(v) = l_Ω(v) + l_Γ(v)

    op = AffineFEOperator(a, l, U, V)
    uh = solve(op)
    writevtk(Λ,"plots/jumps",cellfields=["jump_u"=>jump(uh)])
    writevtk(Ω,"plots/Omega",cellfields=["uh"=>uh])

end


# main2()
main()
