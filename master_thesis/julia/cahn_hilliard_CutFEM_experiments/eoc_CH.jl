##
using Gridap
using Gridap.Algebra
using GridapEmbedded
using Plots
using DataFrames
using CSV
using YAML
using LaTeXStrings
using Latexify
using PrettyTables

function run_CH(;domain="circle", n=2^6, τ_hat, dirname)

    ## Cahn-hilliard
    ε = 1/30
    # ε = 1
    # Gibb's potential
    f(u) = mean(u)*(1 - mean(u)*mean(u))
    u_ex(x, t::Real) = (x[1]*x[1] + x[2]*x[2] - 1 )^4*cos(x[1])*cos(x[2])*cos(t*ε^2)
    u_ex(t) = x -> u_ex(x, t)
    g_0(t) = x -> ( ∂t(u_ex)(x,t) +ε*Δ(Δ(u_ex(t)))(x)
                   # - ( 3/ε )*(2*∇(u_ex(t))(x)⋅∇(u_ex(t))(x) + u_ex(t)(x)*u_ex(t)(x)*Δ(u_ex(t))(x)  )
                  )

    L=2.70
    # n = 2^6
    h = 2*L/n
    # it = 10
    γ = 20
    τ = τ_hat*( ε^2/10 )
    T = ε^2
    it = T/τ
    γg1 = 10
    γg2 = 0.5

    pmin = Point(-L/2, -L/2)
    pmax = Point(L/2, L/2)
    partition = (n,n)
    bgmodel = CartesianDiscreteModel(pmin, pmax, partition)


    graphicsdir = dirname*"/graphics"
    mkpath(graphicsdir)


    # Implicit geometry
    if domain=="circle"
        R  = 1.0
        geo = disk(R)
    elseif domain=="flower"
        function ls_flower(x)
            r0, r1 = L*0.3, L*0.1
            theta = atan(x[1], x[2])
            r0 + r1*cos(5.0*theta) -(x[1]^2 + x[2]^2)^0.5
        end
        # using ! operator to define the interioir
        geo = !AnalyticalGeometry(x-> ls_flower(x))
    end


    # Cut the background model
    cutgeo = cut(bgmodel, geo)
    cutgeo_facets = cut_facets(bgmodel,geo)

    # Set up interpolation mesh and function spaces
    Ω_act = Triangulation(cutgeo, ACTIVE)
    Ω = Triangulation(cutgeo, PHYSICAL)
    Ω_bg = Triangulation(bgmodel)
    writevtk(Ω_bg,   graphicsdir*"/Omega_bg")
    writevtk(Ω,         graphicsdir*"/Omega")
    writevtk(Ω_act,     graphicsdir*"/Omega_act")

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

    function jump_nn(u,n)
        return ( n.plus⋅ ∇∇(u).plus⋅ n.plus - n.minus ⋅ ∇∇(u).minus ⋅ n.minus )
    end

    a_CIP(u,v) = ( ∫(Δ(v)⊙Δ(u))dΩ
                  + ∫(-mean(Δ(v))⊙jump(∇(u)⋅n_Λ) - mean(Δ(u))⊙jump(∇(v)⋅n_Λ) + (γ/h)⋅jump(∇(u)⋅n_Λ)⊙jump(∇(v)⋅n_Λ))dΛ
                  + ∫(-Δ(v)⊙∇(u)⋅n_Γ - Δ(u)⊙∇(v)⋅n_Γ + (γ/h)⋅ ∇(u)⊙n_Γ⋅∇(v)⊙n_Γ )dΓ
                 )

    g(u,v) = h^(-2)*( ∫( (γg1*h)*jump(n_Fg⋅∇(u))*jump(n_Fg⋅∇(v)) ) * dFg +
                     ∫( (γg2*h^3)*jump_nn(u,n_Fg)*jump_nn(v,n_Fg) ) * dFg)

    A_h(u,v) = a_CIP(u,v) + g(u,v)
    lhs(u,v) = ∫(u*v)*dΩ + τ*ε*A_h(u,v)

    c_h(u,v) = ( ∫(f(u)*Δ(v))*dΩ - ∫(f(mean(u))*jump(∇(v)⋅n_Λ))*dΛ - ∫(f(u)*∇(v)⋅n_Γ )*dΓ)
    l_h(v, t ) = ∫(g_0(t)*v)*dΩ
    rhs(u, v, t) =  ∫(u*v)*dΩ + τ*l_h(v,t)# + ( τ/ε) *c_h(u,v)
    rhs(u, t ) = v -> rhs(u,v,t)

    ## time loop
    t0 = 0.0
    Nt_max = convert(Int64, ceil((T - t0)/τ))
    Nt = 0
    t = t0

    # Maximal number of Picard iterations
    kmax = 1

    # Initial data
    u_dof_vals = (rand(Float64, num_free_dofs(U)) .-0.5)*2.0
    uh = interpolate_everywhere(u_ex(0),U)
    # u0 = FEFunction(U, deepcopy(u_dof_vals))
    # uh = FEFunction(U, u_dof_vals)
    # u0_L1 = sum( ∫(u0)dΩ )
    pvd = Dict()
    pvd[t] = createvtk(Ω, graphicsdir*"/sol-tau-$τ_hat-n-$n-$t"*".vtu",cellfields=["uh"=>uh, "u_ex"=>u_ex(t)])

    # Adding initial plotting values
    el2_ts = Float64[]
    eh1_ts = Float64[]
    eah_ts = Float64[]
    ts = Float64[]

    ## Set up linear algebra system
    A = assemble_matrix(lhs, U, V)
    lu = LUSolver()
    cache = nothing

    # Time loop
    println("Solving Cahn-Hilliard for step $(Nt_max), n = $n, τ_hat = $τ_hat")
    while t < T
        Nt += 1
        t += τ
        # println("----------------------------------------")
        k = 0

        if Nt % 25 == 0
            println("$(Nt)/$(Nt_max)")
        end

        while k < kmax
            k += 1
            b = assemble_vector(rhs(uh, t), V)
            op = AffineOperator(A, b)
            cache = solve!(u_dof_vals, lu, op, cache, isnothing(cache))
            uh = FEFunction(U, u_dof_vals)
        end

        # Adding initial plotting values
        # println("----------------------------------------")
        pvd[t] = createvtk(Ω, graphicsdir*"/sol-tau-$τ_hat-n-$n-$t"*".vtu",cellfields=["uh"=>uh, "u_ex"=>u_ex(t)])

        e = u_ex(t) - uh
        el2_t = sqrt(sum( ∫(e*e)dΩ ))
        eh1_t = sqrt(sum( ∫( e*e + ∇(e)⋅∇(e) )*dΩ ))
        eah_t = sqrt(sum( ∫( Δ(e)⋅Δ(e) )*dΩ ))
        push!( ts, t)
        push!( el2_ts, el2_t )
        push!( eh1_ts, eh1_t )
        push!( eah_ts, eah_t )
    end

    # Construct pvd file
    createpvd(graphicsdir*"/sol-tau-$τ_hat-n-$n") do pvd_file
        for (t, vtk) in pvd
            pvd_file[t] = vtk
        end
    end

    parameters = Dict(
        "domain" => domain,
        "gamma" => γ,
        "epsilon" => ε,
        "gamma1" => γg1,
        "gamma2" => γg2,
        "epsilon" => ε,
        "tau_hat" => τ_hat,
        "L"=>2.70,
        "n"=> 2^7,
        "it"=> it
    )

    YAML.write_file(dirname*"/sol-tau-$τ_hat-n-$n.yml", parameters)

    el2_L2 = sqrt(sum(τ* e_ti^2 for e_ti in el2_ts))
    eh1_L2 = sqrt(sum(τ* e_ti^2 for e_ti in eh1_ts))
    eah_L2 = sqrt(sum(τ* e_ti^2 for e_ti in eah_ts))
    el2_inf = maximum(abs.(el2_ts))
    eh1_inf = maximum(abs.(eh1_ts))
    eah_inf = maximum(abs.(eah_ts))
    return el2_L2, eh1_L2, eah_L2, el2_inf, eah_L2, eah_inf

end

function generate_plot(;ns, τs,
        el2s_L2, eh1s_L2, eahs_L2,
        el2s_inf, eh1s_inf, eahs_inf,
        dirname::String, spatial=false, transient=false)

    function compute_eoc(Xs::Vector, errs::Vector)
        eoc = log.(errs[1:end-1] ./ errs[2:end]) ./ log.(Xs[1:end-1] ./ Xs[2:end])
        return [NaN; eoc]
    end

    if spatial
        endfix="spatial"
        Xs = 1 ./ ns
    elseif transient
        endfix="transient"
        Xs = τs
    elseif spatial && transient
        endfix="diagonal"
        Xs = τs
    end

    eoc_el2s_L2 = compute_eoc(Xs, el2s_L2)
    eoc_eh1s_L2 = compute_eoc(Xs, eh1s_L2)
    eoc_eahs_L2 = compute_eoc(Xs, eahs_L2)
    eoc_el2s_inf = compute_eoc(Xs, el2s_inf)
    eoc_eh1s_inf = compute_eoc(Xs, eh1s_inf)
    eoc_eahs_inf = compute_eoc(Xs, eahs_inf)

    # Producing a CSV file
    data = DataFrame( ns=ns, taus=τs,
                     # L2 data
                     el2s_L2=eoc_el2s_L2,
                     eoc_el2s_L2=eoc_el2s_L2,
                     el2s_inf=eoc_el2s_inf,
                     eoc_el2s_inf=eoc_el2s_inf,
                     # H1 data
                     eh1s_L2=eh1s_L2,
                     eoc_eh1s_L2=eoc_eh1s_L2,
                     eh1s_inf=eh1s_inf,
                     eoc_eh1s_inf=eoc_eh1s_inf,
                     # Ah data
                     eahs_L2=eahs_L2,
                     eoc_eahs_L2=eoc_eahs_L2,
                     eahs_inf=eahs_inf,
                     eoc_eahs_inf=eoc_eahs_inf,
                    )
    println(data)
    CSV.write("$dirname/$endfix.csv", delim=',')

end

function conv_test()

    domain="circle"
    maindir = "figures/eoc_CH_$domain"
    if isdir(maindir)
        rm(maindir; recursive=true)
        mkpath(maindir)
    end

    el2s_L2 = Float64[]
    eh1s_L2 = Float64[]
    eahs_L2 = Float64[]
    el2s_inf = Float64[]
    eh1s_inf = Float64[]
    eahs_inf = Float64[]

    # Spatial EOC
    println("Run spatial EOC tests")
    dirname = maindir*"/conv_spatial"
    mkpath(dirname)

    # ns = [2^3, 2^4, 2^5, 2^6, 2^7, 2^8]
    ns = [2^3, 2^4, 2^5, 2^6]
    τ_spatial = 2^(-1)
    for n in ns
        el2_L2, eh1_L2, eah_L2, el2_inf, eh1_inf, eah_inf = run_CH(n=n, τ_hat=τ_spatial, dirname=dirname)
        push!(el2s_L2, el2_L2)
        push!(eh1s_L2, eh1_L2)
        push!(eahs_L2, eah_L2)
        push!(el2s_inf, el2_inf)
        push!(eh1s_inf, eh1_inf)
        push!(eahs_inf, eah_inf)
    end

    τs_spatial = τ_spatial*ones(length(ns))
    generate_plot(ns=ns,τs=τs_spatial,
                  el2s_L2=el2s_L2, eh1s_L2=eh1s_L2, eahs_L2=eahs_L2,
                  el2s_inf=el2s_inf, eh1s_inf=eh1s_inf, eahs_inf=eahs_inf,
                  dirname=maindir, spatial=true)

end

conv_test()

