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

function run_CH(;domain="circle", n=2^6, dirname)

    ## Cahn-hilliard
    ε = 1/30
    # ε = 1
    # Gibb's potential
    f(u) = mean(u)*(1 - mean(u)*mean(u))
    u_ex(x, t::Real) = (x[1]*x[1] + x[2]*x[2] - 1 )^2*cos(x[1])*cos(x[2])*exp(-(4*ε^2 + 2)*t)
    u_ex(t) = x -> u_ex(x, t)
    g_0(t) = x -> ( ∂t(u_ex)(x,t) +ε*Δ(Δ(u_ex(t)))(x)
                   # - ( 3/ε )*(2*∇(u_ex(t))(x)⋅∇(u_ex(t))(x) + u_ex(t)(x)*u_ex(t)(x)*Δ(u_ex(t))(x)  )
                  )

    L=2.70
    # n = 2^6
    h = 2*L/n
    it = 10
    γ = 20
    τ = ε^2/30
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
    rhs(u, v, t) =  ∫(u*v)*dΩ + τ*l_h(v,t) #+ ( τ/ε) *c_h(u,v)
    rhs(u, t ) = v -> rhs(u,v,t)

    ## time loop
    t0 = 0.0
    T = it*τ
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
    pvd[t] = createvtk(Ω, graphicsdir*"/sol-tau-$τ-n-$n-$t"*".vtu",cellfields=["uh"=>uh, "u_ex"=>u_ex(t)])

    # Adding initial plotting values
    el2_ts = Float64[]
    eh1_ts = Float64[]
    ts = Float64[]

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
            b = assemble_vector(rhs(uh, t), V)
            op = AffineOperator(A, b)
            cache = solve!(u_dof_vals, lu, op, cache, isnothing(cache))
            uh = FEFunction(U, u_dof_vals)
        end

        # Adding initial plotting values
        println("----------------------------------------")
        pvd[t] = createvtk(Ω, graphicsdir*"/sol-tau-$τ-n-$n-$t"*".vtu",cellfields=["uh"=>uh, "u_ex"=>u_ex(t)])

        e = u_ex(t) - uh
        el2_t = sqrt(sum( ∫(e*e)dΩ ))
        eh1_t = sqrt(sum( ∫( e*e + ∇(e)⋅∇(e) )*dΩ ))
        println("el2_t $el2_t, eh1_t, $eh1_t ")
        push!( ts, t)
        push!( el2_ts, el2_t )
        push!( eh1_ts, eh1_t )
    end

    # Construct pvd file
    createpvd(graphicsdir*"/sol") do pvd_file
        for (t, vtk) in pvd
            pvd_file[t] = vtk
        end
    end

    parameters = Dict(
        "domain" => domain,
        "gamma" => γ,
        "gamma1" => γg1,
        "gamma2" => γg2,
        "epsilon" => ε,
        "tau" => "epsilon^2/10",
        "L"=>2.70,
        "n"=> 2^7,
        "it"=> it
    )

    YAML.write_file(dirname*"/parameters.yml", parameters)

    el2s_L2 = sqrt(sum(τ* e_ti^2 for e_ti in el2_ts))
    eh1s_L2 = sqrt(sum(τ* e_ti^2 for e_ti in eh1_ts))
    el2s_inf = maximum(abs.(el2_ts))
    eh1s_inf = maximum(abs.(eh1_ts))
    return el2s_L2, eh1s_L2, el2s_inf, eh1s_inf

end

function generate_plot(Xs,
        el2s_L2, eh1s_L2,
        el2s_inf, eh1s_inf,
        dirname::String, Xs_name::String)

    # L2 norms
    p = Plots.plot(Xs, el2s_L2, label="L2L2", legend=:bottomright, xscale=:log2, yscale=:log2, minorgrid=true)
    Plots.scatter!(p, Xs, el2s_L2, primary=false)

    Plots.plot!(p, Xs, eh1s_L2, label=L"L2H1")
    Plots.scatter!(p, Xs, eh1s_L2, primary=false)

    # inf norms
    Plots.plot!(p, Xs, el2s_inf, label=L"infL2")
    Plots.scatter!(p, Xs, el2s_inf, primary=false)

    Plots.plot!(p, Xs, eh1s_inf, label=L"infH1")
    Plots.scatter!(p, Xs, eh1s_inf, primary=false)

    # Configs
    Plots.xlabel!(p, "$Xs_name")
    Plots.plot!(p, xscale=:log2, yscale=:log2, minorgrid=true)
    Plots.plot!(p, legendfontsize=12)  # Adjust the value 12 to your desired font size

    # Save the plot as a .png file using the GR backend
    Plots.savefig(p, "$dirname/plot.png")

    compute_eoc(Xs,  errs) = log.(errs[1:end-1]./errs[2:end])./( log.(Xs[1:end-1]./Xs[2:end]) )
    eoc_l2s_L2 = compute_eoc(Xs, el2s_L2)
    eoc_eh1s_L2 = compute_eoc(Xs, eh1s_L2)
    eoc_l2s_inf = compute_eoc(Xs, el2s_inf)
    eoc_eh1s_inf = compute_eoc(Xs, eh1s_inf)

    eoc_l2s_L2 =  [nothing; eoc_l2s_L2]
    eoc_eh1s_L2 =  [nothing; eoc_eh1s_L2]
    eoc_l2s_inf =  [nothing; eoc_l2s_inf]
    eoc_eh1s_inf =  [nothing; eoc_eh1s_inf]


    Xs_str =  latexify.(Xs)
    data = hcat(Xs_str,
                el2s_L2,  eoc_l2s_L2, eh1s_L2, eoc_eh1s_L2,
                el2s_inf,  eoc_l2s_inf, eh1s_inf, eoc_eh1s_inf)

    formatters = (ft_nonothing, ft_printf("%.1E", [2, 4, 6, 8]), ft_printf("%.2f", [3, 5, 7, 9]))

    header = [Xs_name,
              "L2L2", "EOC", "L2H1", "EOC",
              "infL2", "EOC", "infH1", "EOC"]

    pretty_table(data, header=header, formatters =formatters )
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
    el2s_inf = Float64[]
    eh1s_inf = Float64[]

    # Spatial EOC
    println("Run spatial EOC tests")
    dirname = maindir*"/conv_spatial"
    mkpath(dirname)

    ns = [2^3, 2^4, 2^5, 2^6]
    for n in ns
        el2_L2, eh1_L2, el2_inf, eh1_inf = run_CH(n=n, dirname=dirname)
        push!(el2s_L2, el2_L2)
        push!(eh1s_L2, eh1_L2)
        push!(el2s_inf, el2_inf)
        push!(eh1s_inf, eh1_inf)
    end

    hs = 1 .// ns
    generate_plot(hs,
                  el2s_L2, eh1s_L2,
                  el2s_inf, eh1s_inf,
                  dirname, "h")

end

conv_test()

