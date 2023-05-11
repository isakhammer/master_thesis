
module Solver
    using Test
    # using Plots
    using Dates
    using LaTeXStrings
    using Latexify
    using PrettyTables

    using Gridap
    using Gridap.Algebra

    # ε = 1/100
    ε = 1

    function run(τ, n)
        ## Cahn-hilliard
        # Gibb's potential
        # ψ(u) = u*(1-u^2)
        # u_ex(x, t::Real) = cos(x[1])*cos(x[2])*exp(-4*t*ε^2)
        u_ex(x, t::Real) = cos(x[1])*cos(x[2])*sin(4*t)
        u_ex(t) = x -> u_ex(x,t)
        f(t) = x ->  ∂t(u_ex)(x,t) + Δ(Δ(u_ex(t)))(x)

        ##
        L=2π
        # n = 10
        h = L/n
        domain2D = (0, L, 0, L)
        partition2D = (n, n)
        model = CartesianDiscreteModel(domain2D,partition2D) |> simplexify
        resultdir = "cahn-hilliard-results/"
        mkpath(resultdir)
        writevtk(model, joinpath(resultdir,"model"))

        ## Function spaces
        order = 2
        V = TestFESpace(model, ReferenceFE(lagrangian,Float64, order), conformity=:H1)
        U = TrialFESpace(V)
        Ω = Triangulation(model)
        Λ = SkeletonTriangulation(model)
        Γ = BoundaryTriangulation(model)

        ## Forms
        degree = 2*order
        dΩ = Measure(Ω,degree)
        dΓ = Measure(Γ,degree)
        dΛ = Measure(Λ,degree)

        n_Λ = get_normal_vector(Λ)
        n_Γ = get_normal_vector(Γ)

        # M+ dt A
        γ = 1.5*order*( order+1)
        # τ = ε^2/10
        a(u,v) = ∫(u*v)*dΩ + τ*ε^2*( ∫(Δ(v)⊙Δ(u))dΩ
                                    + ∫(-mean(Δ(v))⊙jump(∇(u)⋅n_Λ) - mean(Δ(u))⊙jump(∇(v)⋅n_Λ) + (γ/h)⋅jump(∇(u)⋅n_Λ)⊙jump(∇(v)⋅n_Λ))dΛ
                                    + ∫(-Δ(v)⊙∇(u)⋅n_Γ - Δ(u)⊙∇(v)⋅n_Γ + (γ/h)⋅ ∇(u)⊙n_Γ⋅∇(v)⊙n_Γ )dΓ
                                   )

        # l(u, v) = ∫(τ*f(t)*v + u*v)*dΩ
        l(u, v) = ∫(u*v)*dΩ
        l(u) = v -> l(u,v)

        ts = Float64[]
        el2_ts = Float64[]
        eh1_ts = Float64[]

        createpvd(resultdir*"ch-solution") do pvd

            ## Set up linear algebra system
            A = assemble_matrix(a, U, V)
            lu = LUSolver()

            ## time loop
            t0 = 0.0
            T = 100*τ
            t = t0

            # u_dof_vals = rand(Float64, num_free_dofs(U))
            uh = interpolate_everywhere(u_ex(0),U)
            pvd[t] = createvtk(Ω, resultdir*"ch-solution_$t"*".vtu",cellfields=["uh"=>uh, "u_ex"=>u_ex(t)])

            # Solve for first time step
            b = assemble_vector(l(uh), V)
            op = AffineOperator(A, b)

            u_dof_vals = get_free_dof_values(uh)
            cache = solve!(u_dof_vals, lu, op)
            uh = FEFunction(U, u_dof_vals)

            t += τ
            pvd[t] = createvtk(Ω, resultdir*"ch-solution_$t"*".vtu",cellfields=["uh"=>uh, "u_ex"=>u_ex(t)])

            # Remaining while loop
            while t < T
                b = assemble_vector(l(uh), V)

                u_dof_vals = get_free_dof_values(uh)
                op = AffineOperator(A, b)
                cache = solve!(u_dof_vals, lu, op, cache, false)
                uh = FEFunction(U, u_dof_vals)

                t += τ
                e = uh - u_ex(t)
                el2_t = sqrt(sum( ∫(e*e)dΩ ))
                eh1_t = sqrt(sum( ∫( e*e + ∇(e)⋅∇(e) )*dΩ ))

                # push!( ts, t)
                push!( el2_ts, el2_t )
                push!( eh1_ts, eh1_t )
                pvd[t] = createvtk(Ω, resultdir*"ch-solution_$t"*".vtu",cellfields=["uh"=>uh, "u_ex"=>u_ex(t),"e"=>e ])
            end

        end

        el2s_L2 = sqrt(sum(τ* e_ti^2 for e_ti in el2_ts))
        eh1s_L2 = sqrt(sum(τ* e_ti^2 for e_ti in eh1_ts))

        el2s_inf = maximum(abs.(el2_ts))
        eh1s_inf = maximum(abs.(eh1_ts))

        return el2s_L2, eh1s_L2, el2s_inf, eh1s_inf
    end


    function generate_figures(Xs::Vector, Xs_name::String,
            el2s_L2::Vector, eh1s_L2::Vector,
            el2s_inf::Vector, eh1s_inf::Vector)

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

        formatters = (ft_nonothing, ft_printf("%.1E", [2, 4, 6, 8]),
                      ft_printf("%.2f", [3, 5, 7, 9]))

        header = [Xs_name,
                  "L2L2", "EOC", "L2H1", "EOC",
                  "infL2", "EOC", "infH1", "EOC"]

        pretty_table(data, header=header, formatters=formatters)
    end


    function spatial()
        ns = [2^2, 2^3, 2^4, 2^5, 2^6]
        τ_const = ε^2/10
        el2s_L2 = Float64[]
        eh1s_L2 = Float64[]
        el2s_inf = Float64[]
        eh1s_inf = Float64[]

        println("SPATIAL TEST")
        for n in ns
            el2_L2,eh1_L2,el2_inf, eh1_inf  = run(τ_const, n)
            println((el2_L2,eh1_L2,el2_inf, eh1_inf))
            push!(el2s_L2, el2_L2)
            push!(eh1s_L2, eh1_L2)
            push!(el2s_inf, el2_inf)
            push!(eh1s_inf, eh1_inf)
        end

        hs = 1 .// ns
        hs_str = "hs"
        generate_figures(hs, "h",
                         el2s_L2, eh1s_L2,
                         el2s_inf, eh1s_inf)
    end

    function diagonal()
        ns = [2^2, 2^3, 2^4, 2^5, 2^6]
        tmp = [ 10^1, 10^2, 10^3, 10^4, 10^5]
        τs_rational = 1 .// tmp # Converting to rational
        τs = ε*ε .*τs_rational
        hs = 1 .// ns
        hs_str = "hs"

        if length(τs) != length( ns )
            return "ERROR"
        end

        println("DIAGONAL TEST")

        el2s_L2 = Float64[]
        eh1s_L2 = Float64[]
        el2s_inf = Float64[]
        eh1s_inf = Float64[]

        for i in 1:length(ns)
            τi = τs[i]
            ni = ns[i]
            el2_L2,eh1_L2,el2_inf, eh1_inf  = run(τi, ni)
            push!(el2s_L2, el2_L2)
            push!(eh1s_L2, eh1_L2)
            push!(el2s_inf, el2_inf)
            push!(eh1s_inf, eh1_inf)
        end

        generate_figures(hs, "h",
                         el2s_L2, eh1s_L2,
                         el2s_inf, eh1s_inf)
        generate_figures(τs_rational, "tau",
                         el2s_L2, eh1s_L2,
                         el2s_inf, eh1s_inf)

    end

end

Solver.diagonal()
Solver.spatial()
