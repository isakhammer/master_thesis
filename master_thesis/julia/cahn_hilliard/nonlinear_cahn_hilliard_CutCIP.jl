
using Test
using Gridap
using Plots
using Dates
using LaTeXStrings
using Latexify
using PrettyTables

module Solver
    using Gridap
    using LinearAlgebra
    using PROPACK
    using GridapEmbedded
    using GridapGmsh
    using Parameters
    import Gridap: ∇

    α = 1
    function man_sol(u_ex)
        # Here is the problem!!
        f(t) = x ->  α*∂t(u_ex)(x,t) + Δ(Δ(u_ex(t)))(x) + 10*u_ex(x,t)*u_ex(x,t)
        ∇u_ex(t) = x ->  ∇(u_ex(t))(x)
        ∇Δu_ex(t) = x ->  ∇(Δ(u_ex(t)))(x)
        return u_ex, f, ∇u_ex, ∇Δu_ex
    end

    @with_kw struct Config
        exact_sol
        circle
    end

    @with_kw struct Solution
        el2s_L2
        eh1s_L2
        ehs_energy_L2
        el2s_inf
        eh1s_inf
        ehs_energy_inf
    end

    function run(;n, dt, solver_config,  vtkdirname=nothing)
        u_ex, f, ∇u_ex, ∇Δu_ex = solver_config.exact_sol
        order=2

        # Background model
        L = 1.11
        domain = (-L, L, -L, L)
        pmin = Point(-L, -L)
        pmax = Point(L, L)
        partition = (n,n)
        bgmodel = CartesianDiscreteModel(pmin, pmax, partition)

        # Implicit geometry
        R  = 1.0
        geo = disk(R)

        # Cut the background model
        cutgeo = cut(bgmodel, geo)
        cutgeo_facets = cut_facets(bgmodel,geo)

        # Set up interpolation mesh and function spaces
        Ω_act = Triangulation(cutgeo, ACTIVE)

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

        # Set up normal vectors
        n_Γ = get_normal_vector(Γ)
        n_Λ = get_normal_vector(Λ)
        n_Fg = get_normal_vector(Fg)

        # Define weak form
        h = L/n

        function mean_n(u,n)
            return 0.5*( u.plus⋅n.plus + u.minus⋅n.minus )
        end

        function mean_nn(u,n)
            return 0.5*( n.plus⋅ ∇∇(u).plus⋅ n.plus + n.minus ⋅ ∇∇(u).minus ⋅ n.minus )
        end

        function jump_nn(u,n)
            return ( n.plus⋅ ∇∇(u).plus⋅ n.plus - n.minus ⋅ ∇∇(u).minus ⋅ n.minus )
        end

        # Define weak form
        γ = 10

        # Ghost penalty parameter
        γg1 = 10/2
        γg2 = 0.1

        # Inner facets
        a(t,u,v) =( ∫( ∇∇(v)⊙∇∇(u) )dΩ
                 + ∫(-mean_nn(v,n_Λ)⊙jump(∇(u)⋅n_Λ) - mean_nn(u,n_Λ)⊙jump(∇(v)⋅n_Λ))dΛ
                 + ∫(-( n_Γ ⋅ ∇∇(v)⋅ n_Γ )⊙∇(u)⋅n_Γ - ( n_Γ ⋅ ∇∇(u)⋅ n_Γ )⊙∇(v)⋅n_Γ)dΓ
                 + ∫((γ/h)⋅jump(∇(u)⋅n_Λ)⊙jump(∇(v)⋅n_Λ))dΛ + ∫((γ/h)⋅ ∇(u)⊙n_Γ⋅∇(v)⊙n_Γ )dΓ
                )

        # Define linear form
        g_1(t) = ∇u_ex(t)⋅n_Γ
        g_2(t) = ∇Δu_ex(t)⋅n_Γ
        b(t,v) = (∫( f(t)*v ) * dΩ
                  +  ∫(-(g_2(t)⋅v))dΓ
                  + ∫(g_1(t)⊙(-(n_Γ⋅∇∇(v)⋅n_Γ) + (γ/h)*∇(v)⋅n_Γ)) * dΓ
               )

        # # Define of ghost penalties
        if order != 2
            println("Not supported order:", order)
        end

        g(t,u,v) = h^(-2)*( ∫( (γg1*h)*jump(n_Fg⋅∇(u))*jump(n_Fg⋅∇(v)) ) * dFg +
                         ∫( (γg2*h^3)*jump_nn(u,n_Fg)*jump_nn(v,n_Fg) ) * dFg)

        F(t,v,u) = ∫( 10*u*u*v )dΩ
        A(t,u,v) = a(t,u,v) + g(t,u,v) + F(t,v,u)
        # Initializing linear terms
        m(t, u, v) = ∫( α* u⋅v )dΩ

        # Alternative 1 (linear)
        op1 = TransientAffineFEOperator(m,A,b,U,V)

        # Alternative 2 (general form)
        res(t,u,v) = A(t,u,v) + m(t, ∂t(u),v) - b(t,v)
        jac(t,u,du,v) = A(t,du,v)
        jac_t(t,u,dut,v) = m(t, dut,v)
        op2 = TransientFEOperator(res,jac,jac_t,U,V)

        # Alternative 3 (Algorithmic differential)
        op3 = TransientFEOperator(res, U, V) # Does not work

        # Intalization of method
        op = op2

        # Solving time problem
        # solver_method = LUSolver()
        solver_method = NLSolver(LUSolver();show_trace=true,method=:newton) #line

        th = 1
        ode_solver = ThetaMethod(solver_method,dt,th)

        s = LUSolver()

        # Inital condition
        t_0 = 0
        T = 1.0
        U_0 = interpolate_everywhere(0,U(0.0))

        #################

        U_h_t = solve(ode_solver, op, U_0, t_0, T)

        ts = Float64[]
        el2_ts = Float64[]
        eh1_ts = Float64[]
        eh_energy_ts = Float64[]

        solname = vtkdirname*"/sol_dt_$dt"*"_n_$n"
        mkpath(solname)
        println("\ndt = ", string(dt), ", n = "*string(n))
        createpvd(solname*".pvd") do pvd
            for (U_h, t) in U_h_t
                e = u_ex(t) - U_h
                pvd[t] = createvtk(Ω, solname*"_$t"*".vtu",cellfields=["u_h"=>U_h,"e"=>e])
                el2_t = sqrt(sum( ∫(e*e)dΩ ))
                eh1_t = sqrt(sum( ∫( e*e + ∇(e)⋅∇(e) )*dΩ ))

                eh_energy = sqrt(sum( ∫(e⊙e)*dΩ + ∫( ∇∇(e)⊙∇∇(e) )*dΩ
                              + ( γ/h ) * ∫(jump(∇(e)⋅n_Λ) ⊙ jump(∇(e)⋅n_Λ))dΛ
                              + ( h/γ ) * ∫(mean_nn(e,n_Λ) ⊙ mean_nn(e,n_Λ))dΛ
                              + ( γ/h ) * ∫((∇(e)⋅n_Γ) ⊙ (∇(e)⋅n_Γ))dΓ
                              + ( h/γ ) * ∫(( n_Γ ⋅ ∇∇(e)⋅ n_Γ ) ⊙ ( n_Γ ⋅ ∇∇(e)⋅ n_Γ ))dΓ
                             ))
                push!( ts, t)
                push!( el2_ts, el2_t )
                push!( eh1_ts, eh1_t )
                push!( eh_energy_ts, eh_energy)
            end
        end

        el2s_L2 = sqrt(sum(dt* e_ti^2 for e_ti in el2_ts))
        eh1s_L2 = sqrt(sum(dt* e_ti^2 for e_ti in eh1_ts))
        ehs_energy_L2 = sqrt(sum(dt* e_ti^2 for e_ti in eh_energy_ts))

        el2s_inf = maximum(abs.(el2_ts))
        eh1s_inf = maximum(abs.(eh1_ts))
        ehs_energy_inf = maximum(abs.(eh_energy_ts))

        sol = Solution(el2s_L2=el2s_L2, eh1s_L2=eh1s_L2, ehs_energy_L2=ehs_energy_L2,
                       el2s_inf=el2s_inf, eh1s_inf=eh1s_inf, ehs_energy_inf=ehs_energy_inf)
        return sol
    end

end # Solver

function generate_figures(Xs,
        el2s_L2, eh1s_L2, ehs_energy_L2,
        el2s_inf, eh1s_inf, ehs_energy_inf,
        dirname::String, Xs_name::String)

    filename = dirname*"/conv_"*Xs_name
    compute_eoc(Xs, errs) = log.(errs[1:end-1]./errs[2:end])./log.(Xs[1:end-1]./Xs[2:end])


    Xs_str =  latexify.(Xs)

    eoc_l2s_L2 = compute_eoc(Xs, el2s_L2)
    eoc_eh1s_L2 = compute_eoc(Xs, eh1s_L2)
    eoc_ehs_energy_L2 = compute_eoc(Xs, ehs_energy_L2)
    eoc_l2s_inf = compute_eoc(Xs, el2s_inf)
    eoc_eh1s_inf = compute_eoc(Xs, eh1s_inf)
    eoc_ehs_energy_inf = compute_eoc(Xs, ehs_energy_inf)

    eoc_l2s_L2 =  [nothing; eoc_l2s_L2]
    eoc_eh1s_L2 =  [nothing; eoc_eh1s_L2]
    eoc_ehs_energy_L2 =  [nothing; eoc_ehs_energy_L2]
    eoc_l2s_inf =  [nothing; eoc_l2s_inf]
    eoc_eh1s_inf =  [nothing; eoc_eh1s_inf]
    eoc_ehs_energy_inf =  [nothing; eoc_ehs_energy_inf]

    data = hcat(Xs_str,
                el2s_L2,  eoc_l2s_L2, eh1s_L2, eoc_eh1s_L2, ehs_energy_L2, eoc_ehs_energy_L2,
                el2s_inf,  eoc_l2s_inf, eh1s_inf, eoc_eh1s_inf, ehs_energy_inf, eoc_ehs_energy_inf)

    formatters = ( ft_nonothing, ft_printf("%.2f", [3, 5, 7, 9, 11,13]),
                  ft_printf("%.1E", [2, 4, 6, 8, 10,12]))

    header = ["$Xs_name",
              "L2L2", "EOC", "L2H1", "EOC", "L2ah", "EOC",
              "infL2", "EOC", "infH1", "EOC", "infah", "EOC"]

    pretty_table(data, header=header, formatters =formatters )


    formatters = ( ft_nonothing, ft_printf("%.5f", [3, 5, 7, 9, 11,13]),
                  ft_printf("%.1E", [2, 4, 6, 8, 10,12]))

    open(filename*".tex", "w") do io
        pretty_table(io, data, header=header, backend=Val(:latex ), formatters = formatters )
    end

    # L2 norms
    p = Plots.plot(Xs, el2s_L2, label="L2L2", legend=:bottomright, xscale=:log2, yscale=:log2, minorgrid=true)
    Plots.scatter!(p, Xs, el2s_L2, primary=false)

    Plots.plot!(p, Xs, eh1s_L2, label=L"L2H1")
    Plots.scatter!(p, Xs, eh1s_L2, primary=false)

    Plots.plot!(p, Xs, ehs_energy_L2, label=L"L2ah")
    Plots.scatter!(p, Xs, ehs_energy_L2, primary=false)

    # inf norms
    Plots.plot!(p, Xs, el2s_inf, label=L"infL2")
    Plots.scatter!(p, Xs, el2s_inf, primary=false)

    Plots.plot!(p, Xs, eh1s_inf, label=L"infH1")
    Plots.scatter!(p, Xs, eh1s_inf, primary=false)

    Plots.plot!(p, Xs, ehs_energy_inf, label=L"infah")
    Plots.scatter!(p, Xs, ehs_energy_inf, primary=false)

    # Configs
    Plots.xlabel!(p, "$Xs_name")
    Plots.plot!(p, xscale=:log2, yscale=:log2, minorgrid=true)
    Plots.plot!(p, legendfontsize=12)  # Adjust the value 12 to your desired font size

    # Save the plot as a .png file using the GR backend
    # Plots.gr()
    # Plots.pgfplotsx()
    Plots.savefig(p, filename*"_plot.png")
    # Plots.savefig(p, filename*"_plot.tex")

end
function convergence_analysis(; ns, dts, dirname, solver_config, spatial=false, dt_const=2^-4, transient=false, n_const=2^7)
    println("Run convergence",)

    if (transient)
        # Time dim EOC
        el2s_L2 = Float64[]
        eh1s_L2 = Float64[]
        ehs_energy_L2 = Float64[]
        el2s_inf = Float64[]
        eh1s_inf = Float64[]
        ehs_energy_inf = Float64[]
        println("Run convergence tests with constant n = "*string(n_const))

        for dt in dts
            sol = Solver.run(n=n_const, dt=dt, solver_config=solver_config, vtkdirname=dirname)

            push!(el2s_L2, sol.el2s_L2)
            push!(eh1s_L2, sol.eh1s_L2)
            push!(ehs_energy_L2, sol.ehs_energy_L2)
            push!(el2s_inf, sol.el2s_inf)
            push!(eh1s_inf, sol.eh1s_inf)
            push!(ehs_energy_inf, sol.ehs_energy_inf)
        end
        generate_figures(dts, el2s_L2, eh1s_L2, ehs_energy_L2,
                         el2s_inf, eh1s_inf, ehs_energy_inf,
                         dirname, "dt")
    end


    if (spatial)
        # Spatial EOC
        el2s_L2 = Float64[]
        eh1s_L2 = Float64[]
        ehs_energy_L2 = Float64[]
        el2s_inf = Float64[]
        eh1s_inf = Float64[]
        ehs_energy_inf = Float64[]
        println("Run convergence tests with constant dt = "*string(dt_const))
        for n in ns
            sol = Solver.run(n=n, dt=dt_const, solver_config=solver_config, vtkdirname=dirname)
            push!(el2s_L2, sol.el2s_L2)
            push!(eh1s_L2, sol.eh1s_L2)
            push!(ehs_energy_L2, sol.ehs_energy_L2)
            push!(el2s_inf, sol.el2s_inf)
            push!(eh1s_inf, sol.eh1s_inf)
            push!(ehs_energy_inf, sol.ehs_energy_inf)
        end
        hs =  1 .// ns
        generate_figures(hs, el2s_L2, eh1s_L2, ehs_energy_L2,
                         el2s_inf, eh1s_inf, ehs_energy_inf,
                         dirname, "h")
    end



end

function main_convergence()

    dirname= "figures/nonlinear_cahn_hilliard_CutCIP/conv_"*string(Dates.now())
    println(dirname)
    mkpath(dirname)

    u_ex(x,t::Real) = sin(t)*(x[1]^2 + x[2]^2 - 1 )^3*sin(x[1])*cos(x[2])
    u_ex(t::Real) = x -> u_ex(x,t)
    exact_sol = Solver.man_sol(u_ex)
    circle = true
    solver_config = Solver.Config(exact_sol, circle)

    dts = [2^-2,2^-3,2^-4,2^-5]
    ns = [2^2,2^3,2^4,2^5, 2^6]

    @time convergence_analysis( ns=ns, dts=dts, dirname=dirname, solver_config=solver_config, transient=false, spatial=true)



end

function main_simulation_test()
    dirname= "figures/nonlinear_cahn_hilliard_CutCIP/sim_"*string(Dates.now())
    println(dirname)
    mkpath(dirname)

    u_ex(x,t::Real) = sin(t)*(x[1]^2 + x[2]^2 - 1 )^3*sin(x[1])*cos(x[2])
    u_ex(t::Real) = x -> u_ex(x,t)
    exact_sol = Solver.man_sol(u_ex)
    circle = true
    solver_config = Solver.Config(exact_sol, circle)
    @time Solver.run(n=2^5, dt=2^-3, solver_config=solver_config, vtkdirname=dirname)
end

main_simulation_test()

