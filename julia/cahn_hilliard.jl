
using Test
using Gridap
using Plots
using LaTeXStrings
# plotlyjs()


function run_cahn_hilliard(; n=10, order::Int, generate_vtk=false, dirname="cahn_hilliard", test=false)

    # Analytical manufactured solution
    α = 1
    u(x) = cos(x[1])*cos(x[2])
    f(x) = Δ(Δ(u))(x)+ α*u(x)
    g(x) = Δ(u)(x)


    # Domain
    L = 2*π

    @test f(VectorValue(0.5,0.5)) == ( 4+α )*u(VectorValue(0.5,0.5))
    @test g(VectorValue(0.5,0.5)) == -2*u(VectorValue(0.5,0.5))

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

    b(t, v) = ∫( v*f )dΩ + ∫( g*(∇(v)⋅n_Γ) )dΓ

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
    U_h_t = solve(ode_solver, op, U_0, t_0, T)

    solname = dirname*"/sol"
    createpvd(solname) do pvd
        for (U_h, t) in U_h_t
            println("t "*string(t))
            pvd[t] = createvtk(Ω, solname*"_$t"*".vtu",cellfields=["U_h"=>U_h])
        end
    end



    # Error

    # e = u - uh
    # l2(u) = sqrt(sum( ∫( u⊙u )*dΩ ))
    # h1(u) = sqrt(sum( ∫( u⊙u + ∇(u)⊙∇(u) )*dΩ ))
    # el2 = l2(e)
    # eh1 = h1(e)

    # if test==true
    #     @test el2 < 10^-3
    # end

    # if !generate_vtk
    #     return el2, eh1
    # end

    # writevtk(model, dirname*"/model")
    # writevtk(Λ,dirname*"/skeleton")
    # writevtk(Λ,dirname*"/jumps",cellfields=["jump_u"=>jump(uh)])
    # writevtk(Ω,dirname*"/omega",cellfields=["uh"=>uh])
    # writevtk(Ω,dirname*"/error",cellfields=["e"=>e])
    # writevtk(Ω,dirname*"/manufatured",cellfields=["u"=>u])

    # return el2, eh1

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
