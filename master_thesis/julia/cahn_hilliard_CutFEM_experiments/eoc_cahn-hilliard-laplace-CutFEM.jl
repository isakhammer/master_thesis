##
using Gridap
using Gridap.Algebra
using GridapEmbedded

function main()
    ## Cahn-hilliard
    ε = 1/30
    f(u) = mean(u)*(mean(u)*mean(u) - 1) # why mean necessary???
    # f(u) = u*(u*u - 1)
    f_der(u) = 3*u*u - 1

    u_ex(x, t::Real) = cos(x[1])*cos(x[2])*exp(-(4*ε^2 + 2)*t)
    u_ex(t) = x -> u_ex(x,t)

    ##
    L=2.50
    n = 2^7
    h = L/n
    pmin = Point(-L, -L)
    pmax = Point(L, L)
    partition = (n,n)
    bgmodel = CartesianDiscreteModel(pmin, pmax, partition)

    resultdir = "CutFEM-cahn-hilliard-results/"
    mkpath(resultdir)
    writevtk(bgmodel, joinpath(resultdir,"model"))

    # Implicit geometry
    domain = "other"
    if domain=="circle"
        R  = 1.0
        geo1 = disk(R)
    elseif domain=="other"
        geo1 = disk(0.5, x0=Point(0.0,0.0))
        geo2 = disk(0.5, x0=Point(0.95,0.0))
        geo3 = disk(0.5, x0=Point(0.95,0.95))
        geo4 = disk(0.5, x0=Point(0.0,0.95))

        geo = union(geo1,geo2)
        geo = union(geo,geo3)
        geo = union(geo,geo4)
    end


    # Cut the background model
    cutgeo = cut(bgmodel, geo)
    cutgeo_facets = cut_facets(bgmodel,geo)

    # Set up interpolation mesh and function spaces
    Ω_act = Triangulation(cutgeo, ACTIVE)

    ## Function spaces
    order = 2
    # Construct function spaces
    V = TestFESpace(Ω_act, ReferenceFE(lagrangian, Float64, order), conformity=:H1)
    U = TrialFESpace(V)


    # Set up integration meshes, measures and normals
    Ω = Triangulation(cutgeo, PHYSICAL)
    Γ = EmbeddedBoundary(cutgeo)
    Λ = SkeletonTriangulation(cutgeo_facets)
    Fg = GhostSkeleton(cutgeo)

    # Set up integration measures
    degree = 2*order
    dΩ   = Measure(Ω, degree)
    dΓ   = Measure(Γ, degree)
    dΛ   = Measure(Λ, degree) # F_int
    dFg  = Measure(Fg, degree)

    n_Λ = get_normal_vector(Λ)
    n_Γ = get_normal_vector(Γ)

    # Set up normal vectors
    n_Γ = get_normal_vector(Γ)
    n_Λ = get_normal_vector(Λ)
    n_Fg = get_normal_vector(Fg)

    # Manufactured solution
    g_0(t) = x -> ( ∂t(u_ex)(x,t) +ε*Δ(Δ(u_ex(t)))(x)
                   - ( 3/ε )*(2*∇(u_ex(t))(x)⋅∇(u_ex(t))(x) + u_ex(t)(x)*u_ex(t)(x)*Δ(u_ex(t))(x)  )
                  )

    ∇u_ex(t) = x ->  ∇(u_ex(t))(x)
    g_1(t) = ∇u_ex(t)⋅n_Γ

    ∇Δu_ex(t) = x ->  ∇(Δ(u_ex(t)))(x)
    g_2(t) = ∇Δu_ex(t)⋅n_Γ


    function jump_nn(u,n)
        return ( n.plus⋅ ∇∇(u).plus⋅ n.plus - n.minus ⋅ ∇∇(u).minus ⋅ n.minus )
    end

    # M+ dt A
    γ = 1.5*order*( order+1)
    τ = ε^2/10
    a_L(u,v) =   ( ∫(Δ(v)⊙Δ(u))dΩ
                    + ∫(-mean(Δ(v))⊙jump(∇(u)⋅n_Λ) - mean(Δ(u))⊙jump(∇(v)⋅n_Λ) + (γ/h)⋅jump(∇(u)⋅n_Λ)⊙jump(∇(v)⋅n_Λ))dΛ
                    + ∫(-Δ(v)⊙∇(u)⋅n_Γ - Δ(u)⊙∇(v)⋅n_Γ + (γ/h)⋅ ∇(u)⊙n_Γ⋅∇(v)⊙n_Γ )dΓ
                   )

    γg1 = 10/2
    γg2 = 0.1
    g_h(u,v) = h^(-2)*( ∫( (γg1*h)*jump(n_Fg⋅∇(u))*jump(n_Fg⋅∇(v)) ) * dFg +
                     ∫( (γg2*h^3)*jump_nn(u,n_Fg)*jump_nn(v,n_Fg) ) * dFg)

    a(u,v) = a_L(u,v) + g_h(u,v)

    c_h(u,v) = ( ∫(-f(u)*Δ(v))*dΩ
                + ∫(f(u)*jump(∇(v)⋅n_Λ))*dΛ
                + ∫(f(u)*∇(v)⋅n_Γ )*dΓ
                + ∫(f_der(u)*g_1(t)*∇(v)⋅n_Γ )*dΓ
               )

    l_L(t,v) = ∫(g_0(t)*v)*dΩ - ∫(ε*g_2(t)* v)*dΓ - ∫(ε*g_1(t)*Δ(v))*dΓ

    # Constructing  right hand side  (which is explicit)
    rhs(t, u, v) = τ*l_L(t,v) + ∫(u*v)*dΩ + τ *(1/ε)*c_h(u,v)
    rhs(t, u) = v -> rhs(t,u,v)

    lhs(u,v) = ∫(u*v)*dΩ  + τ*ε*a(u,v)
    lhs(u,v) =  τ*ε*a(u,v)


    ## time loop
    t0 = 0.0
    N_it = 100
    T = N_it*τ
    Nt_max = convert(Int64, ceil((T - t0)/τ))
    Nt = 0
    t = t0

    # Maximal number of Picard iterations
    kmax = 1

    # Initial data
    u_dof_vals = (rand(Float64, num_free_dofs(U)) .-0.5)*2.0
    # uh = interpolate_everywhere(u_ex(0),U)
    uh = FEFunction(U, u_dof_vals)
    pvd = Dict()
    pvd[t] = createvtk(Ω, resultdir*"ch-solution_$t"*".vtu",cellfields=["uh"=>uh, "u_ex"=>u_ex(t)])

    println("========================================")
    println("Solving Cahn-Hilliard with t0 = $t0, T = $T and time step τ = $τ with Nt_max = $Nt_max timesteps")
    println("========================================")

    ## Set up linear algebra system
    A = assemble_matrix(lhs, U, V)
    lu = LUSolver()
    cache = nothing

    # Time loop
    while t < T
        Nt += 1
        t += τ
        println("----------------------------------------")
        println("Solving Cahn-Hilliard for t = $t, step $(Nt)/$(Nt_max)")
        k = 0
        while k < kmax
            k += 1
            println("Iteration k = $k")
            b = assemble_vector(rhs(t, uh), V)
            op = AffineOperator(A, b)
            cache = solve!(u_dof_vals, lu, op, cache, isnothing(cache))
            uh = FEFunction(U, u_dof_vals)
        end
        println("----------------------------------------")
        pvd[t] = createvtk(Ω, resultdir*"ch-solution_$t"*".vtu",cellfields=["uh"=>uh, "u_ex"=>u_ex(t)])
    end

    # Construct pvd file
    createpvd(resultdir*"ch-solution") do pvd_file
        for (t, vtk) in pvd
            pvd_file[t] = vtk
        end
    end



end

main()

