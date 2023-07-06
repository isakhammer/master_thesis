##
using Gridap
using Gridap.Algebra
using Gridap.CellData
using Gridap.Visualization
using GridapEmbedded
using Plots
using DataFrames
using CSV
using YAML
using LaTeXStrings
using Random
using Dates
Random.seed!(1234)  # Set the seed to a specific value

# function main(;ε::Float64 = 1/100, domain="flower", L::Float64=1.1, τ=1/10^5, n::Int64=2^6, it::Int64=10, result_dir_suffix::String ="")
function main(;ε = 1/100, domain="flower", L=2.70, τ=1/10^5, n=2^6, it=10, result_dir_suffix ="")

    ## Cahn-hilliard
    # ε = 1/100
    # ε = 1
    # Gibb's potential
    f(u) = mean(u)*(mean(u)*mean(u) - 1)

    ##
    # L= 2.70
    # n = 2^6
    h = 2*L/n
    # it = 10000
    γ = 20
    # τ = ε^2/60
    γg1 = 10
    γg2 = 0.5

    pmin = Point(-L/2, -L/2)
    pmax = Point(L/2, L/2)
    partition = (n,n)
    bgmodel = CartesianDiscreteModel(pmin, pmax, partition)

    # Define maindir, add result_dir_suffix if not trivial
    maindir = "figures"*(s-> all(isspace, s) ? "" : "_$(s)")(result_dir_suffix)*"/physical_CH_$(domain)"
    if isdir(maindir)
        rm(maindir; recursive=true)
        mkpath(maindir)
    end

    graphicsdir = maindir*"/graphics"
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

    # Refined model for plotting P2 properly

    partition_ref = (2n,2n)
    bgmodel_ref = CartesianDiscreteModel(pmin, pmax, partition_ref)
    cutgeo_ref = cut(bgmodel_ref, geo)
    # Set up interpolation mesh and function spaces
    Ω_act_ref = Triangulation(cutgeo_ref, ACTIVE)
    Ω_ref = Triangulation(cutgeo_ref, PHYSICAL)
    Ω_bg_ref = Triangulation(bgmodel_ref)
    writevtk(Ω_bg_ref,   graphicsdir*"/Omega_bg_ref")
    writevtk(Ω_ref,         graphicsdir*"/Omega_ref")
    writevtk(Ω_act_ref,     graphicsdir*"/Omega_act_ref")

    ## Function spaces
    order = 2
    # Construct function spaces
    V = TestFESpace(Ω_act, ReferenceFE(lagrangian, Float64, order), conformity=:H1)
    U = TrialFESpace(V)

    # Construct refined spaces for output
    V_ref = FESpace(Ω_act_ref, ReferenceFE(lagrangian, Float64, order), conformity=:H1)

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
    lhs(u,v) = ∫(u*v)*dΩ + τ*ε^2*A_h(u,v)
    # Properly scaled ghost penalty including mass matrix scaling
    # lhs(u,v) = ∫(u*v)*dΩ + τ*ε^2*a_CIP(u,v) + (τ*ε^2+h^4)*g(u,v)

    c_h(u,v) = ( ∫(f(u)*(-1)*Δ(v))*dΩ + ∫(f(mean(u))*jump(∇(v)⋅n_Λ))*dΛ + ∫(f(u)*∇(v)⋅n_Γ )*dΓ)
    rhs(u, v) =  ∫(u*v)*dΩ -τ*c_h(u,v)
    rhs(u) = v -> rhs(u,v)

    # Defining plotting values
    δuhs = []
    Δuhs = []
    Es = Float64[]
    E1s = Float64[]
    E2s = Float64[]
    ts = Float64[]

    createpvd(graphicsdir*"/sol") do pvd

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
        u0 = FEFunction(U, u_dof_vals)

        # u_ex(x) = sin(2π*x[1])*sin(2π*x[2])
        # u0 = interpolate_everywhere(u_ex, U)

        # Initializing discrete function and its dof cals
        uh = FEFunction(U,deepcopy( u0.free_values ) )
        iuh = Interpolable(uh)
        uh_ref = interpolate_everywhere(iuh, V_ref)

        # defining
        u0_L1 = abs( sum( ∫(u0)dΩ ) )
        # pvd[t] = createvtk(Ω, graphicsdir*"/sol_$t"*".vtu",cellfields=["uh"=>uh])
        pvd[t] = createvtk(Ω_ref, graphicsdir*"/sol_$t"*".vtu",cellfields=["uh"=>uh_ref])

        # Adding initial plotting values
        E = sum( ∫(( ε^2/2 )*( ∇(uh)⋅∇(uh) ) + (1/4)*((uh*uh - 1)*(uh*uh - 1))  )dΩ)
        E1 = sum( ∫(( ε^2/2*∇(uh)⋅∇(uh) ) )dΩ)
        E2 = sum( ∫((1/4)*((uh*uh - 1)*(uh*uh - 1))  )dΩ)
        push!(ts, t)
        push!(Es, E)
        push!(E1s, E1)
        push!(E2s, E2)
        # First step is not defined
        push!( δuhs, missing)
        push!( Δuhs, missing)


        println("========================================")
        println("Solving Cahn-Hilliard with ε = $ε, t0 = $t0, T = $T and time step τ = $τ with Nt_max = $Nt_max timesteps")
        println("========================================")

        ## Set up linear algebra system
        A = assemble_matrix(lhs, U, V)
        lu = LUSolver()
        cache = nothing

        # Time loop
        println("----------------------------------------")
        println("Solving Cahn-Hilliard for $domain with step $(Nt_max), n = $n, τ = $τ")

        uh0 = FEFunction(U, deepcopy( uh.free_values ))
        while t < T
            Nt += 1
            t += τ
            k = 0
            if Nt%10 == 0
                println("$(Nt)/$(Nt_max)")
            end

            while k < kmax
                k += 1
                b = assemble_vector(rhs(uh), V)
                op = AffineOperator(A, b)
                # Updating the dof values to uh
                cache = solve!(uh.free_values, lu, op, cache, isnothing(cache))
            end

            # Adding plotting values
            E = sum( ∫(ε^2/2*( ∇(uh)⋅∇(uh) ) + (1/4)*((uh*uh - 1)*(uh*uh - 1))  )dΩ)
            E1 = sum( ∫((ε^2/2*∇(uh)⋅∇(uh) ) )dΩ)
            E2 = sum( ∫((1/4)*((uh*uh - 1)*(uh*uh - 1))  )dΩ)
            push!(ts, t)
            push!(Es, E)
            push!(E1s, E1)
            push!(E2s, E2)
            δuh =  sum( ∫( (uh0 - uh ) )dΩ) /u0_L1
            Δuh =  abs(sum( ∫( (u0 - uh ) )dΩ)) /u0_L1
            push!( δuhs, δuh)
            push!( Δuhs, Δuh)

            # Updating previous step
            uh0 = FEFunction(U, deepcopy( uh.free_values ))
            iuh = Interpolable(uh)
            uh_ref = interpolate_everywhere(iuh, V_ref)
            # pvd[t] = createvtk(Ω, graphicsdir*"/sol_$t"*".vtu",cellfields=["uh"=>uh])
            pvd[t] = createvtk(Ω_ref, graphicsdir*"/sol_$t"*".vtu",cellfields=["uh"=>uh_ref])
        end
    end

    # Save results
    df = DataFrame(ts=ts, Es=Es, E1s=E1s, E2s=E2s, delta_uhs=δuhs, Delta_uhs=Δuhs)
    CSV.write(maindir*"/sol.csv", df, delim=',')

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
    YAML.write_file(maindir*"/parameters.yml", parameters)

    its = LinRange(1, length(ts), length(ts))

    default_size = (800, 1000)
    p1 = plot(its[1:end-1], δuhs[2:end], size=default_size, xscale=:log10, legend=false, ylabel=L"$\delta u$", xlabel=L"$t/\tau$")
    scatter!(its[1:end-1], δuhs[2:end], markersize = 2)
    p2 = plot(its[1:end-1], Δuhs[2:end], size=default_size, yscale=:log10, legend=false, xscale=:log10, ylabel=L"$\Delta u$", xlabel=L"$t/\tau$")
    scatter!(its[1:end-1], Δuhs[2:end], markersize = 2)
    p3 = plot(its, Es, size=default_size, yscale=:log10, xscale=:log10, legend=false, ylabel=L"$E^m$", xlabel=L"$t/\tau$")
    scatter!(its, Es, markersize = 2)
    p4 = plot(its[1:end-1], Es[1:end-1] .- Es[2:end], size=default_size, legend=false,  xscale=:log10,  xlabel=L"$t/\tau$", ylabel=L"$\Delta E^m$")
    scatter!(its[1:end-1], Es[1:end-1] .- Es[2:end], markersize = 2)
    p5 = plot(its, E1s, size=default_size, yscale=:log10, xscale=:log10, legend=false, ylabel=L"$E_1^m$", xlabel=L"$t/\tau$")
    scatter!(its, E1s, markersize = 2)
    p6 = plot(its, E2s, size=default_size, xscale=:log10, legend=false, ylabel=L"$E_2^m$", xlabel=L"$t/\tau$")
    scatter!(its, E2s, markersize = 2)
    p = plot(p1, p2, p3, p4, p5, p6, layout = (6, 1))

    savefig(p,maindir*"/physical_CH_plot.pdf" )
    plot(p)

end

ε = 1/100
domain="circle"
L=2.70
τ = ε^2
n=2^6
it=10000
# result_dir_suffix = string(Dates.now())
result_dir_suffix = ""
# result_dir_suffix = "test"
main(;ε, domain="circle", L, τ, n, it, result_dir_suffix)
main(;ε, domain="flower", L, τ, n, it, result_dir_suffix)
# function main(;ε = 1/100, domain="flower", L=1.1, τ=1/10^5, n=2^6, it=10, result_dir_suffix ="")
# main(;domain)
