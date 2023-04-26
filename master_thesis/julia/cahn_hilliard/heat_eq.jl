

using Test
using Plots
using Dates
using LaTeXStrings
using Latexify
using PrettyTables

module HE
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

    function run(;n::Number, dt::Number,
            vtkdirname::String, ode_method::String="BE", u_ex::Union{Function, Nothing}=nothing)

        u_ex, ∇u_ex, f = man_sol(u_ex)
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

        # Iterative solvers
        algebraic_solver = LUSolver()
        # solver_method = NLSolver(LUSolver();show_trace=true,method=:newton) #line

        # ODE solvers

        function ode_solver(ode_method, algebraic_solver, dt)
            if ode_method == "BE"
                return ThetaMethod(algebraic_solver, dt, 1)
            elseif ode_method == "CN"
                return ThetaMethod(algebraic_solver, dt, 0.5)
            else
                error("Invalid solver_choice: $ode_method")
            end
        end

        solver = ode_solver(ode_method, algebraic_solver, dt)

        # Inital condition
        t_0 = 0
        T = 1.0
        U_0 = interpolate_everywhere(0, U(0.0))

        #################

        U_h_t = solve(solver, op, U_0, t_0, T)

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


function main()
    vtkdirname= "figures/heat_eq/test_"*string(Dates.now())
    println(vtkdirname)
    mkpath(vtkdirname)
    u_ex(x,t::Real) = sin(t)*(x[1]^2 + x[2]^2 - 1 )^3*sin(x[1])*cos(x[2])
    u_ex(t::Real) = x -> u_ex(x,t)
    n, dt = 2^6, 2^(-4)
    HE.run(n=n, dt=dt, vtkdirname=vtkdirname, ode_method="BE", u_ex=u_ex)
end

# main_simulation_test()
# main()
