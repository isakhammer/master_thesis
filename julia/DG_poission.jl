using Gridap

function main2()
    # source: https://gridap.github.io/Tutorials/v0.5/pages/t005_dg_discretization/
    u(x) = 3*x[1] + x[2]^2
    f(x) = -2 - 12*x[1]
    g(x) = u(x)

    ∇u(x) = VectorValue(3   ,
                        2*x[2] )

    # mesh generation
    L = 1.0
    n = 4
    order = 2

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

    # Generate plots
    dirname = "possion_results"
    if (isdir(dirname))
        rm(dirname, recursive=true)
    end
    mkdir(dirname)

    writevtk(Λ, dirname*"/poission_skeleton")
    writevtk(Λ, dirname*"/poission_jumps", cellfields=["jump_u"=>jump(uh)])
    writevtk(Γ, dirname*"/poission_shell", cellfields=["uh"=>uh])
    writevtk(Ω, dirname*"/poission_field", cellfields=["uh"=>uh])

end

main2()

