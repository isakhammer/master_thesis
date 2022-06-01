using Gridap
import Gridap: ∇


function main2()
    # source: https://gridap.github.io/Tutorials/v0.5/pages/t005_dg_discretization/
    u(x) = 3*x[1] + x[2]^2
    f(x) = -2 - 12*x[1]
    g(x) = u(x)

    # ∇u(x) = VectorValue(3   + x[2]*x[3],
    #                 2*x[2]   + x[1]*x[3],
    #                 6*x[3]^2 + x[1]*x[2])


    println(∇(u) === ∇u)

    # mesh generation
    L = 1.0
    n = 4
    order = 2

    domain2D = (0.0, L, 0.0, L)
    partition2D = (n,n)
    model = CartesianDiscreteModel(domain2D,partition2D)
    writevtk(model,"model_mesh")

    V = TestFESpace(model, ReferenceFE(lagrangian,Float64,order), conformity=:L2)
    U = TrialFESpace(V)

    Ω = Triangulation(model)
    Γ = BoundaryTriangulation(model)
    Λ = SkeletonTriangulation(model)

    writevtk(Λ,"strian_cube")

    degree = 2*order
    dΩ = Measure(Ω,degree)
    dΓ = Measure(Γ,degree)
    dΛ = Measure(Λ,degree)

    n_Γ = get_normal_vector(Γ)
    n_Λ = get_normal_vector(Λ)

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
    writevtk(Λ, "jumps", cellfields=["jump_u"=>jump(uh)])
    writevtk(Γ, "shell", cellfields=["uh"=>uh])
    writevtk(Ω, "field", cellfields=["uh"=>uh])

end

function main()
    f(x) = 1
    g(x) = u(x)

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
    writevtk(model,"model_mesh")

    V = TestFESpace(model, ReferenceFE(lagrangian,Float64,order), conformity=:L2)
    U = TrialFESpace(V)
    Ω = Triangulation(model)
    Γ = BoundaryTriangulation(model)

    Λ = SkeletonTriangulation(model)

    writevtk(Λ,"strian_cube")

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
    writevtk(Λ,"jumps",cellfields=["jump_u"=>jump(uh)])
    writevtk(Ω,"Omega",cellfields=["uh"=>uh])

end


main2()
# main()
