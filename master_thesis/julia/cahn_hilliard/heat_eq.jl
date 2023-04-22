

using Test
using Plots
using Dates
using LaTeXStrings
using Latexify
using PrettyTables

module Solver
    using Gridap
    using LinearAlgebra
    using PROPACK
    using Gridap.ODEs.ODETools: GeneralizedAlpha
    using Gridap.ODEs.ODETools: Newmark
    using GridapEmbedded
    using GridapGmsh
    using Parameters
    import Gridap: ∇

    function man_sol(u_ex)
        α = 1
        f(t) = x ->  ∂t(u_ex)(x,t) - Δ(u_ex(t))(x)
        ∇u_ex(t) = x ->  ∇(u_ex(t))(x)
        return u_ex, ∇u_ex, f
    end

    @with_kw struct Config
        exact_sol
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
        u_ex, ∇u_ex, f = solver_config.exact_sol
        order = 2

        L= 1.0
        domain = (0,L,0,L)
        partition = (n,n)
        model = CartesianDiscreteModel(domain, partition)
        reffe = ReferenceFE(lagrangian, Float64, order)
        V = FESpace( model, reffe, conformity=:H1)
        U = TrialFESpace(V)

        Ω = Triangulation(model)
        Γ = BoundaryTriangulation(model)
        degree = 2*order
        dΩ = Measure(Ω,degree)
        dΓ = Measure(Γ,degree)
        n_Γ = get_normal_vector(Γ)

        g(t) = ∇u_ex(t)⋅n_Γ
        a(u,v) = ∫(∇(v)⋅∇(u))dΩ
        b(v,t) = ∫(v*f(t))dΩ + ∫((g(t)⋅v))dΓ

        res(t,u,v) = a(u,v) + ∫(∂t(u)*v)dΩ - b(v,t)
        jac(t,u,du,v) = a(du,v)
        jac_t(t,u,dut,v) = ∫(dut*v)dΩ

        op = TransientFEOperator(res,jac,jac_t,U,V)

        # Nonlinear/Linear solver
        solver_method =LUSolver()
        # solver_method = NLSolver(LUSolver();show_trace=true,method=:newton) #linesearch=BackTracking())

        # ODE solver
        ode_solver = ThetaMethod(solver_method,dt, 1) # Works!
        # ode_solver = RungeKutta(solver_method,dt,:BE_1_0_1) # Works!
        # ode_solver = RungeKutta(solver_method,dt,:SDIRK_2_1_2) # Does not converge!
        # ode_solver = RungeKutta(solver_method,dt,:TRBDF2_3_3_2) # Does not converge!
        # γ, β  = 0.5, 0.25
        # ode_solver = Newmark(solver_method,dt,γ,β) # Does not compile
        # ρ∞ = 1.0 # Equivalent to Newmark(0.5, 0.25)
        # ode_solver = GeneralizedAlpha(solver_method, dt, ρ∞) # Does not compile

        # Inital condition
        t_0 = 0
        T = 1.0
        U_0 = interpolate_everywhere(0, U(0.0))

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
                eh_energy = sqrt(sum( ∫(∇(e)⋅∇(e) )*dΩ ))
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

    dirname= "figures/heat_eq/conv_"*string(Dates.now())
    println(dirname)
    mkpath(dirname)

    u_ex(x,t) = sin(t)*cos(2π*x[1])*sin(2π*x[2])
    u_ex(t::Real) = x -> u_ex(x,t)
    exact_sol = Solver.man_sol(u_ex)
    solver_config = Solver.Config(exact_sol)

    dts = [2^-2,2^-3,2^-4,2^-5]
    ns = [2^2,2^3,2^4,2^5, 2^6]
    @time convergence_analysis( ns=ns, dts=dts, dirname=dirname, solver_config=solver_config, transient=true, spatial=true)

end

function main_simulation_test()
    dirname= "figures/heat_eq/sim_"*string(Dates.now())
    println(dirname)
    mkpath(dirname)

    u_ex(x,t) = sin(t)*cos(2π*x[1])*sin(2π*x[2])
    u_ex(t::Real) = x -> u_ex(x,t)
    exact_sol = Solver.man_sol(u_ex)
    solver_config = Solver.Config(exact_sol)
    @time Solver.run(n=2^5, dt=2^-4, solver_config=solver_config, vtkdirname=dirname)
end


main_convergence()
# main_simulation_test()
