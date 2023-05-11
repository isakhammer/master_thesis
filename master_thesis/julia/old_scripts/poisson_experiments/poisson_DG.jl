include("results.jl")
using Dates


module Solver
    using Gridap
    using Parameters

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


    function run(; order=order, n=n, use_quads=true)

        L = 2π
        pmin = Point(0.,0.)
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
        u(x) = 100*cos(x[1])*cos(x[2])
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
        return sol
    end

end # module



function convergence_analysis(;orders, ns, dirname)
    println("Run convergence",)

    for order in orders

        el2s = Float64[]
        eh1s = Float64[]
        ehs_energy = Float64[]
        println("Run convergence tests: order = "*string(order))

        for n in ns

            res = Solver.run(order=order, n=n)

            push!(el2s, res.el2)
            push!(eh1s, res.eh1)
            push!(ehs_energy, res.eh_energy)
        end
        Results.generate_figures(ns=ns, el2s=el2s, eh1s=eh1s, ehs_energy=ehs_energy, order=order, dirname=dirname)
    end

end

function main()
    dirname = "figures/poisson_DG/"*string(Dates.now())
    mkpath(dirname)
    orders=[1,2,3,4]
    ns = [2^2, 2^3, 2^4, 2^5, 2^6, 2^7]
    convergence_analysis(orders=orders, ns=ns, dirname=dirname)
end

@time main()

