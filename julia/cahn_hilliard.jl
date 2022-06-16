
using Test
using Gridap
using Plots
using LaTeXStrings
# plotlyjs()


function run_cahn_hilliard(; n=10, order::Int, generate_vtk=false, dirname="cahn_hilliard", test=false)

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
    Δt = 0.1
    th = 0.5
    ode_solver = ThetaMethod(linear_solver,Δt,th)

    # Inital condition
    U_0 = interpolate_everywhere(0,U(0.0))

    #################


    op = op_Af
    # op = op_AD
    t_0=0
    T=3
    tol = 10^-3

    U_h_t = solve(ode_solver, op, U_0, t_0, T)

    solname = dirname*"/sol"
    createpvd(solname) do pvd
        for (U_h, t) in U_h_t
            println("t "*string(t))
            pvd[t] = createvtk(Ω, solname*"_$t"*".vtu",cellfields=["U_h"=>U_h])

            e = u(t) - U_h
            l2(w) = w*w
            el2 = sqrt(sum( ∫(l2(e))dΩ ))
            @test el2 < tol
        end
    end

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

main()
