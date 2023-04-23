using Dates
using Gridap

function run(;n, u_ex,  dt, solver_choice="ThetaMethod", vtkdirname=nothing)
    # Manufactured solution
    f(t) = x ->  ∂t(u_ex)(x,t) - Δ(u_ex(t))(x)
    ∇u_ex(t) = x ->  ∇(u_ex(t))(x)


    # Setup
    L= 1.0
    domain = (0,L,0,L)
    partition = (n,n)
    model = CartesianDiscreteModel(domain, partition)
    order = 2
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
    s =LUSolver()
    # s = NLSolver(LUSolver();show_trace=true,method=:newton) #linesearch=BackTracking())

    function choose_solver(solver_choice, solver_method, dt)
        if solver_choice == "ThetaMethod"
            return ThetaMethod(solver_method, dt, 1)
        elseif solver_choice == "RungeKutta_BE_1_0_1"
            return RungeKutta(solver_method, dt, :BE_1_0_1)
        elseif solver_choice == "RungeKutta_SDIRK_2_1_2"
            return RungeKutta(solver_method, dt, :SDIRK_2_1_2)
        elseif solver_choice == "RungeKutta_TRBDF2_3_3_2"
            return RungeKutta(solver_method, dt, :TRBDF2_3_3_2)
        elseif solver_choice == "Newmark"
            γ, β  = 0.5, 0.25
            return Newmark(solver_method,dt,γ,β)
        elseif solver_choice == "GeneralizedAlpha"
            ρ∞ = 1.0 # Equivalent to Newmark(0.5, 0.25)
            return GeneralizedAlpha(solver_method, dt, ρ∞) # Does not compile
        else
            error("Invalid solver_choice: $solver_choice")
        end
    end

    ode_solver = choose_solver(solver_choice, s, dt)

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

    println("\nTesting solver choice $solver_choice")
    println("L2 Norms:")
    println("el2s_L2: $(round(el2s_L2, digits=4))")
    println("eh1s_L2: $(round(eh1s_L2, digits=4))")
    println("ehs_energy_L2: $(round(ehs_energy_L2, digits=4))")

    println("Infinity Norms:")
    println("el2s_inf: $(round(el2s_inf, digits=4))")
    println("eh1s_inf: $(round(eh1s_inf, digits=4))")
    println("ehs_energy_inf: $(round(ehs_energy_inf, digits=4))")
end


function main()
    dirname= "figures/heat_eq/sim_"*string(Dates.now())
    println(dirname)
    mkpath(dirname)

    u_ex(x,t) = sin(t)*cos(2π*x[1])*sin(2π*x[2])
    u_ex(t::Real) = x -> u_ex(x,t)
    # run(n=2^5, dt=2^-4, u_ex=u_ex, vtkdirname=dirname, solver_choice="ThetaMethod")
    # run(n=2^5, dt=2^-4, u_ex=u_ex, vtkdirname=dirname, solver_choice="RungeKutta_BE_1_0_1")
    # run(n=2^5, dt=2^-4, u_ex=u_ex, vtkdirname=dirname, solver_choice="RungeKutta_SDIRK_2_1_2")
    # run(n=2^5, dt=2^-4, u_ex=u_ex, vtkdirname=dirname, solver_choice="RungeKutta_TRBDF2_3_3_2")
    run(n=2^5, dt=2^-4, u_ex=u_ex, vtkdirname=dirname, solver_choice="Newmark")
    run(n=2^5, dt=2^-4, u_ex=u_ex, vtkdirname=dirname, solver_choice="GeneralizedAlpha")
end


main()
