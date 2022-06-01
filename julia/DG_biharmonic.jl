using Gridap
import Gridap: ∇

function main()

    u(x) = 3*x[1] + x[2]^2 + 2*x[3]^3 + x[1]*x[2]*x[3]
    f(x) = -2 - 12*x[3]
    g(x) = u(x)

    ∇u(x) = VectorValue(3        + x[2]*x[3],
                        2*x[2]   + x[1]*x[3],
                        6*x[3]^2 + x[1]*x[2])

    L = 1.0
    domain = (0.0, L, 0.0, L)
    n = 4
    partition = (n,n)
    model = CartesianDiscreteModel(domain,partition)

    order = 3
    V = TestFESpace(model,
                    ReferenceFE(lagrangian,Float64,order),
                    conformity=:L2)

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


    writevtk(Λ,"SkeletonTriangulation")
    writevtk(Γ,"BoundaryTriangulation")
    writevtk(Ω,"Triangulation")
    writevtk(model,"model")


    a_Ω(u,v) = ∫( ∇(v)⊙∇(u) )dΩ
    l_Ω(v) = ∫( v*f )dΩ

    h = L / n
    γ = order*(order+1)
    a_Γ(u,v) = ∫( - v*(∇(u)⋅n_Γ) - (∇(v)⋅n_Γ)*u + (γ/h)*v*u )dΓ
    l_Γ(v)   = ∫(                - (∇(v)⋅n_Γ)*g + (γ/h)*v*g )dΓ

    a_Λ(u,v) = ∫( - jump(v*n_Λ)⊙mean(∇(u))
                  - mean(∇(v))⊙jump(u*n_Λ)
                  + (γ/h)*jump(v*n_Λ)⊙jump(u*n_Λ) )dΛ


    a(u,v) = a_Ω(u,v) + a_Γ(u,v) + a_Λ(u,v)
    l(v) = l_Ω(v) + l_Γ(v)

    op = AffineFEOperator(a, l, U, V)
    uh = solve(op)

    writevtk(Ω,"results",cellfields=["uh"=>uh])

end

function main2()
    u(x) = 3*x[1] + x[2]^2 + 2*x[3]^3 + x[1]*x[2]*x[3]
    f(x) = -2 - 12*x[3]
    g(x) = u(x)

    ∇u(x) = VectorValue(3   + x[2]*x[3],
                    2*x[2]   + x[1]*x[3],
                    6*x[3]^2 + x[1]*x[2])


    println(∇(u) === ∇u)

    # mesh generation
    L = 1.0
    domain = (0.0, L, 0.0, L, 0.0, L)
    n = 4
    partition = (n,n,n)
    model = CartesianDiscreteModel(domain,partition)
    writevtk(model,"model_mesh")

    domain2D = (0.0, L, 0.0, L)
    partition2D = (n,n)
    model2D = CartesianDiscreteModel(domain2D,partition2D)

    order = 3
    V = TestFESpace(model,
                ReferenceFE(lagrangian,Float64,order),
                conformity=:L2)
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
    a_Γ(u,v) = ∫( - v*(∇(u)⋅n_Γ) - (∇(v)⋅n_Γ)*u + (γ/h)*v*u )dΓ
    l_Γ(v)   = ∫(                - (∇(v)⋅n_Γ)*g + (γ/h)*v*g )dΓ

    a_Λ(u,v) = ∫( - jump(v*n_Λ)⊙mean(∇(u))
              - mean(∇(v))⊙jump(u*n_Λ)
              + (γ/h)*jump(v*n_Λ)⊙jump(u*n_Λ) )dΛ


    a(u,v) = a_Ω(u,v) + a_Γ(u,v) + a_Λ(u,v)
    l(v) = l_Ω(v) + l_Γ(v)

    op = AffineFEOperator(a, l, U, V)
    uh = solve(op)
    writevtk(Λ,"jumps",cellfields=["jump_u"=>jump(uh)])
    writevtk(Ω,"field",cellfields=["uh"=>uh])





end

main2()
