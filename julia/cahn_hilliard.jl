
using Test
using Gridap
using Plots
using LaTeXStrings
# plotlyjs()


function run_cahn_hilliard(; n=10, order::Int, generate_vtk=false, dirname="cahn_hilliard", test=false, Δt=0.1, t_0=0, T=3)

    # Analytical manufactured solution
    α = 1
    u(x,t::Real) = t*cos(x[1])*cos(x[2])
    u(t::Real) = x -> u(x,t)

    # f(t::Real) = x -> Δ(Δ(u(t)))(x)+ α*∂t(u)
    # g(t::Real) = x -> Δ(u(t))(x)

    # f(x, t::Real) =  1
    # g(x, t::Real) =  0
    # f(t) = x -> ∂t(u)(x,t)-Δ(u(t))(x)

    f(t) = x -> α*∂t(u)(x,t)+Δ(Δ(u(t)))(x)
    g(t) = x ->  Δ(u(t))(x)

    # Manufactured solution test
    # tₜ, x1ₜ, x2ₜ  = 0.5, 0.5, 0.5
    # @test f(tₜ)(VectorValue(x1ₜ,x2ₜ)) ==  4*tₜ*cos(x1ₜ)*cos(x2ₜ) + cos(x1ₜ)*cos(x2ₜ)
    # @test g(tₜ)(VectorValue(x1ₜ,x2ₜ)) == -2*tₜ*cos(x1ₜ)*cos(x2ₜ)

    # u(x,t) = (1.0-x[1])*x[1]*(1.0-x[2])*x[2]*t
    # u(t::Real) = x -> u(x,t)
    # f(t) = x -> ∂t(u)(x,t)-Δ(u(t))(x)


    # u(x) = cos(x[1])*cos(x[2])
    # f(x) = Δ(Δ(u))(x)+ α*u(x)
    # g(x) = Δ(u)(x)


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
    # op = op_AD
    tol = 10^-3

    U_h_t = solve(ode_solver, op, U_0, t_0, T)

    # solname = dirname*"/sol"
    # createpvd(solname) do pvd
    #     for (U_h, t) in U_h_t
    #         println("t "*string(t))
    #         pvd[t] = createvtk(Ω, solname*"_$t"*".vtu",cellfields=["U_h"=>U_h])

    #         e = u(t) - U_h
    #         l2(w) = w*w
    #         el2 = sqrt(sum( ∫(l2(e))dΩ ))
    #         @test el2 < tol
    #     end
    # end

    return model, u, U_h_t, Ω, dΩ

end

function construct_pvd_mansol(dirname, u, U_h_t, Ω)

    pvddirname =dirname*"/pvd"
    if (!isdir(pvddirname))
        mkdir(pvddirname)
    end

    solname = pvddirname*"/sol"

    # mansolname = pvddirname*"/mansol"
    # createpvd(mansolname) do pvd
    #     for (U_h, t) in U_h_t
    #         pvd[t] = createvtk(Ω, mansolname*"_$t"*".vtu",cellfields=["u"=>u(t)])
    #     end
    # end

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

    # Generate plots
    function makedir(dirname)
        if (isdir(dirname))
            rm(dirname, recursive=true)
        end
        mkdir(dirname)
    end

    folder = "cahn_hilliard"
    makedir(folder)

    exampledir = folder*"/example"
    makedir(exampledir)
    run_cahn_hilliard(n=90, order=2, generate_vtk=true, dirname=exampledir, test=false)
end


function searchsortednearest(a,x)
    # https://discourse.julialang.org/t/findnearest-function/4143/5
    idx = searchsortedfirst(a,x)
    if (idx==1); return idx; end
    if (idx>length(a)); return length(a); end
    if (a[idx]==x); return idx; end
    if (abs(a[idx]-x) < abs(a[idx-1]-x))
        return idx
    else
        return idx-1
    end
end


function unit_square_convergence()
    dirname = "cahn_hilliard/unit_square_convergence"

    if (isdir(dirname))
        rm(dirname, recursive=true)
    end
    mkdir( dirname )

    Δt, t_0, T = 0.1, 0, 3

    ns = [ 32, 64, 90]
    # ns = [ 8, 16]#, 32, 64, 90]

    N = size(ns)[1]
    hs = zeros(N)
    ts = Float64[]
    el2_matrix = zeros( ( N, Int(T/Δt) ) )
    eh1_matrix = zeros( ( N, Int(T/Δt) ) )

    for i in 1:N
        n = ns[i]
        println("runs n", n)
        n_dirname = dirname*"/man_sol_se_"*string(n)
        model, u, U_h_t, Ω, dΩ = run_cahn_hilliard(; n=n,  order=2, generate_vtk=false, dirname="cahn_hilliard", test=false, Δt=Δt, t_0=t_0, T=T)
        ts2, el2_ts, eh1_ts = analyze(n_dirname, model, u, U_h_t, Ω, dΩ)
        ts =ts2
        el2_matrix[i,:] = el2_ts
        eh1_matrix[i,:] = eh1_ts
        hs[i] = 2*π/n
    end

    t1 = t_0
    j1 = searchsortednearest(ts,t1)

    t2 = 1.0
    j2 = searchsortednearest(ts,t2)

    t3 = T
    j3 = searchsortednearest(ts,t3)


    Plots.plot(hs, [ el2_matrix[:, j1] el2_matrix[:, j2] el2_matrix[:, j3] ] ,
        yaxis=:log,
        xaxis=:log,
        label=["t = "*string(t1) "t = "*string(t2) "t = "*string(t3)  ],
        shape=:auto,
        xlabel="h",
        ylabel="L2 - error norm")

    Plots.png(dirname*"/L2_convergence_plots")

    Plots.plot(hs, [ eh1_matrix[:, j1] eh1_matrix[:, j2] eh1_matrix[:, j3] ] ,
        yaxis=:log,
        xaxis=:log,
        label=["t = "*string(t1) "t = "*string(t2) "t = "*string(t3)  ],
        shape=:auto,
        xlabel="h",
        ylabel="H1 - error norm")

    Plots.png(dirname*"/H1_convergence_plots")

    println("Here")
    return


end

unit_square_convergence()
# main()
