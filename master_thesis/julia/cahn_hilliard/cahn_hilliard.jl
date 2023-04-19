
using Test
using Gridap
using Plots
using Dates
using LaTeXStrings
# plotlyjs()

# Analytical manufactured solution
α = 1
u(x,t::Real) = t*cos(x[1])*cos(x[2])
u(t::Real) = x -> u(x,t)


f(t) = x -> α*∂t(u)(x,t)+Δ(Δ(u(t)))(x)
g(t) = x ->  Δ(u(t))(x)

# Manufactured solution test
tₜ, x1ₜ, x2ₜ  = 0.5, 0.5, 0.5
@test f(tₜ)(VectorValue(x1ₜ,x2ₜ)) ==  4*tₜ*cos(x1ₜ)*cos(x2ₜ) + cos(x1ₜ)*cos(x2ₜ)
@test g(tₜ)(VectorValue(x1ₜ,x2ₜ)) == -2*tₜ*cos(x1ₜ)*cos(x2ₜ)


function run_cahn_hilliard(; n=10, order::Int, generate_vtk=false, dirname="cahn_hilliard", test=false, Δt=0.1, t_0=0, T=3)

    # Domain
    L = 2*π
    domain = (0,L,0,L)
    partition = (n,n)
    model = CartesianDiscreteModel(domain,partition)

    # FE space
    reffe = ReferenceFE(lagrangian,Float64,order)
    V = TestFESpace(model,reffe)
    U = TransientTrialFESpace(V)

    # Triangulation
    Ω = Triangulation(model)
    Γ = BoundaryTriangulation(model)
    Λ = SkeletonTriangulation(model)
    degree = 2*order
    dΩ = Measure(Ω,degree)
    dΓ = Measure(Γ,degree)
    dΛ = Measure(Λ,degree)

    n_Γ = get_normal_vector(Γ)
    n_Λ = get_normal_vector(Λ)

    # Weak form
    h = (domain[2]-domain[1]) / partition[1]
    h = L / n
    γ = 1

    m(t, u, v) = ∫( α* u⋅v )dΩ
    a(t, u, v) = ∫( Δ(u)*Δ(v) )dΩ +
                ∫( - mean(Δ(u))*jump(∇(v)⋅n_Λ) - jump(∇(u)⋅n_Λ)*mean(Δ(v))
                    + γ/h*jump(∇(u)⋅n_Λ)*jump(∇(v)⋅n_Λ) )dΛ

    b(t, v) = ∫( v*f(t) )dΩ + ∫( g(t)*(∇(v)⋅n_Γ) )dΓ

    # Initializing linear terms
    op_Af = TransientAffineFEOperator(m,a,b,U,V)

    # Solving time problem
    linear_solver = LUSolver()
    th = 0.5
    ode_solver = ThetaMethod(linear_solver,Δt,th)

    # Inital condition
    U_0 = interpolate_everywhere(0,U(0.0))

    #################

    op = op_Af
    U_h_t = solve(ode_solver, op, U_0, t_0, T)

    return model, u, U_h_t, Ω, dΩ

end

function analyze(dirname::String, model, u, U_h_t, Ω, dΩ)

    if (isdir(dirname))
        rm(dirname, recursive=true)
    end
    mkdir(dirname)

    println("Anayzing data and generating plots")

    ts = Float64[]
    el2_ts = Float64[]
    eh1_ts = Float64[]

    solname = dirname*"/sol"
    createpvd(solname) do pvd
        for (U_h, t) in U_h_t
            println("t = "*string(t))
            pvd[t] = createvtk(Ω, solname*"_$t"*".vtu",cellfields=["u_h"=>U_h])
            e = u(t) - U_h
            el2_t = sqrt(sum( ∫(e*e)dΩ ))
            eh1_t = sqrt(sum( ∫( e*e + ∇(e)⋅∇(e) )*dΩ ))

            push!( ts, t)
            push!( el2_ts, el2_t )
            push!( eh1_ts, eh1_t )
        end
    end

    # Checks if all error values are acceptable
    @test all(el2_ts.<10^-2)
    @test all(eh1_ts.<10^-1)


    writevtk(model,dirname*"/model")
    writevtk(Ω,dirname*"/Omega_triangulation", )

    Plots.plot(ts, [ el2_ts eh1_ts ],
        yaxis=:log,
        label=["L2" "H1"],
        shape=:auto,
        xlabel="t",
        ylabel="error norm")

    Plots.png(dirname*"/convergence_plots_se")
    return ts, el2_ts, eh1_ts
end

function main()

    dirname= "figures/CIP_cahn_hilliard/example"*string(Dates.now())
    println(dirname)
    mkpath(dirname)

    model, u, U_h_t, Ω, dΩ = run_cahn_hilliard(n=90, order=2, generate_vtk=true, dirname=dirname, test=false)
    ts, el2_ts, eh1_ts = analyze(dirname, model, u, U_h_t, Ω, dΩ)
end

main()
