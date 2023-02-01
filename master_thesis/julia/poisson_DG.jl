include("results.jl")
using Dates

module Solver
    using Gridap
    using Parameters
    import GridapMakie
    import Makie
    import GLMakie
    using Test

    function man_sol(;L=1,m=1,r=1)
        u(x) = 100*cos(m*( 2π/L )*x[1])*cos(r*( 2π/L )*x[2])
    end


    @with_kw struct Solution
        Ω
        Γ
        Λ

        model
        h::Real

        u
        uh
        e
        el2
        eh1
        eh_energy
    end

    function generate_vtk(; sol::Solution, dirname::String)
        println("Generating vtk's in ", dirname)
        if (isdir(dirname))
            rm(dirname, recursive=true)
        end
        mkdir(dirname)

        writevtk(sol.model, dirname*"/model")
        writevtk(sol.Λ, dirname*"/skeleton")
        writevtk(sol.Γ, dirname*"/boundary")
        writevtk(sol.Λ, dirname*"/jumps",cellfields=["jump_u"=>jump(sol.uh)])
        writevtk(sol.Ω, dirname*"/omega",cellfields=["uh"=>sol.uh])
        writevtk(sol.Ω, dirname*"/error",cellfields=["e"=>sol.e])
        writevtk(sol.Ω, dirname*"/manufatured",cellfields=["u"=>sol.u])


        fig = Makie.plot(sol.Λ)
        Makie.wireframe!(sol.Λ, color=:black, linewidth=2)
        Makie.wireframe!(sol.Γ, color=:black, linewidth=2)
        Makie.save(dirname*"/grid.png", fig)

        # (Isak): Doesnt work :( Please fix
        # fig, _ , plt = Makie.plot(res.Ω, res.uh)
        # Makie.Colorbar(fig[1,2], plt)
        # Makie.save(dirname*"/man_sol.png", fig)
    end


    function run(; order=order, n=n, L=L,  m=m, r=r, dirname=nothing, use_quads=true)

        pmin = Point(0.,0.0)
        pmax = Point(L, L)
        partition = (n, n)

        # h = norm((pmax-pmin)./VectorValue(partition)) # unstable!!!
        h = L / n

        if !use_quads
            model = CartesianDiscreteModel(pmin, pmax, partition) |> simplexify
        else
            model = CartesianDiscreteModel(pmin, pmax, partition)
        end

        # u is the manufactured solution
        # -Δu = f in Ω, and u = g on Γ
        # u(x) = 100*cos(x[1])*cos(x[2])
        u = man_sol(L=L, m=m, r=r)
        # u(x) = 3*x[1] + x[2]^2
        f(x) = -Δ(u)(x)
        g(x) = u(x)

        # Define triangulation
        Ω = Triangulation(model)
        Γ = BoundaryTriangulation(model)
        Λ = SkeletonTriangulation(model)

        ## Function spaces
        reffe = ReferenceFE(lagrangian, Float64, order)

        V = TestFESpace(Ω, reffe, conformity=:L2)
        U = TrialFESpace(V)

        ## Define the weak form
        degree = 2*order
        dΩ = Measure(Ω, degree)
        dΓ= Measure(Γ, degree)
        dΛ= Measure(Λ, degree)

        n_Γ  = get_normal_vector(Γ)
        n_Λ  = get_normal_vector(Λ)

        a_Ω(u,v) =∫( ∇(v)⋅∇(u) )dΩ
        l_Ω(v) = ∫( v⊙f )dΩ

        γ = order*(order+1)  # Penalty parameter
        μ = γ/h

        a_Γ(u,v) =∫( - ( ∇(u)⋅n_Γ )⊙v - u⊙( ∇(v)⋅n_Γ ) + μ*u⊙v )dΓ # isak
        l_Γ(v) = ∫( -(( ∇(v)⋅n_Γ )⊙g) + μ*(g⊙v) )dΓ # isak

        # Comment: Seems like gridap does not like mean(∇(u)⋅n_Λ )
        a_Λ(u,v) =∫( - mean(∇(u))⊙jump(v⋅n_Λ) - jump(u⋅n_Λ)⊙mean(∇(v)) + μ* jump(u)⊙jump(v)  )dΛ

        a(u,v) = a_Ω(u,v) + a_Γ(u,v) + a_Λ(u,v)
        l(v) = l_Ω(v) + l_Γ(v)

        op = AffineFEOperator(a, l, U, V)
        uh = solve(op)

        e = u - uh
        function mean_n(u,n)
            return 0.5*( u.plus⋅n.plus + u.minus⋅n.minus )
        end

        el2 = sqrt(sum( ∫(e*e)dΩ ))
        eh1 = sqrt(sum( ∫( e⊙e + ∇(e)⊙∇(e) )*dΩ ))

        # Defined at eq 2.18 GurkanMassing2019
        eh_energy = sqrt(sum( ∫(∇(e)⊙∇(e) )*dΩ
                             + ∫(h*jump( e )⊙jump( e ) )dΛ
                             # + ∫(h*mean( ∇(e)⋅n_Λ )⊙mean( ∇(e)⋅n_Λ ) )dΛ # does not work :(
                             + ∫(h*mean_n( ∇(e),n_Λ )⊙mean_n( ∇(e),n_Λ ) )dΛ
                             + ∫(h*( ∇(e)⋅n_Γ )⊙( ∇(e)⋅n_Γ ) )dΓ
                            ))

        u_inter = interpolate(u, V)

        sol = Solution(  model=model, Ω=Ω, Γ=Γ, Λ=Λ, h=h,
                        u=u_inter, uh=uh, e=e, el2=el2, eh1=eh1, eh_energy=eh_energy)
        if ( dirname!=nothing)
            generate_vtk(sol, dirname)
        end
        return sol
    end

end # module



function convergence_analysis(; L, m, r, orders, ns, dirname, optimize=false)
    println("Run convergence",)

    for order in orders

        el2s = Float64[]
        eh1s = Float64[]
        ehs_energy = Float64[]
        println("Run convergence tests: order = "*string(order))

        for n in ns

            res = Solver.run(order=order, n=n, L=L, m=m, r=r)

            if !(optimize)
                vtkdirname =dirname*"/order_"*string(order)*"_n_"*string(n)
                mkpath(vtkdirname)
                res = Solver.run(order=order, n=n, L=L, m=m, r=r, dirname=vtkdirname)
            else
                res = Solver.run(order=order, n=n, L=L, m=m, r=r)
            end

            push!(el2s, res.el2)
            push!(eh1s, res.eh1)
            push!(ehs_energy, res.eh_energy)
        end
        Results.generate_figures(ns=ns, el2s=el2s, eh1s=eh1s, ehs_energy=ehs_energy, order=order, dirname=dirname)
    end

end

function main()
    function makedir(dirname)
        if (isdir(dirname))
            rm(dirname, recursive=true)
        end
        mkdir(dirname)
    end

    resultdir= "figures/poisson_DG/"*string(Dates.now())
    mkpath(resultdir)

    function run(;  L,m,r)
        orders=[2,3,4]
        ns = [2^2, 2^3, 2^4, 2^5, 2^6, 2^7]
        dirname = resultdir*"/L_"*string(round(L,digits=2))*"_m_"*string(m)*"_r_"*string(r);
        makedir(dirname)
        convergence_analysis( L=L, m=m, r=r, orders=orders, ns=ns, dirname=dirname, optimize=true)
    end

    run(L=1,m=1,r=1)
end

@time main()

