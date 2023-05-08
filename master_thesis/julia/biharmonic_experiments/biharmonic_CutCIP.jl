using Dates
import Plots

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

    # α(x) = x[1]^2 + x[2]^2
    α(x) = 1

    function man_sol(u_ex)
        ∇u_ex(x) = ∇(u_ex)(x)
        ∇Δu_ex(x) = ∇(Δ(u_ex))(x)
        f(x) = Δ(Δ(u_ex))(x)+ α(x)⋅u_ex(x)
        return u_ex, f, ∇u_ex, ∇Δu_ex
    end

    @with_kw struct Solution
        el2
        eh1
        eh_energy
        cond_number
        ndof
    end


    function run(; order, n, u_ex, dirname)

        println("sim order ", order, " n ", n)
        u_ex, f, ∇u_ex, ∇Δu_ex = man_sol(u_ex)

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
        γ = 10*(order*(order - 1)/2)

        # Ghost penalty parameter
        # γg0 = γ
        γg1 = 10/2
        γg2 = 0.1

        # Mesh size
        h = L/n

        function mean_n(u,n)
            return 0.5*( u.plus⋅n.plus + u.minus⋅n.minus )
        end

        function mean_nn(u,n)
            return 0.5*( n.plus⋅ ∇∇(u).plus⋅ n.plus + n.minus ⋅ ∇∇(u).minus ⋅ n.minus )
        end

        function jump_nn(u,n)
            return ( n.plus⋅ ∇∇(u).plus⋅ n.plus - n.minus ⋅ ∇∇(u).minus ⋅ n.minus )
            # return jump( n⋅ ∇∇(u)⋅ n)
        end

        # Inner facets
        a(u,v) =( ∫( ∇∇(v)⊙∇∇(u) + α⋅(v⊙u) )dΩ
                 + ∫(-mean_nn(v,n_Λ)⊙jump(∇(u)⋅n_Λ) - mean_nn(u,n_Λ)⊙jump(∇(v)⋅n_Λ))dΛ
                 + ∫(-( n_Γ ⋅ ∇∇(v)⋅ n_Γ )⊙∇(u)⋅n_Γ - ( n_Γ ⋅ ∇∇(u)⋅ n_Γ )⊙∇(v)⋅n_Γ)dΓ
                 + ∫((γ/h)⋅jump(∇(u)⋅n_Λ)⊙jump(∇(v)⋅n_Λ))dΛ + ∫((γ/h)⋅ ∇(u)⊙n_Γ⋅∇(v)⊙n_Γ )dΓ
                )

        # Define linear form
        # Notation: g_1 = ∇u_ex⋅n_Γ, g_2 = ∇Δu_ex⋅n_Γ
        g_1 = ∇u_ex⋅n_Γ
        g_2 = ∇Δu_ex⋅n_Γ
        l(v) = (∫( f*v ) * dΩ
                +  ∫(-(g_2⋅v))dΓ
                + ∫(g_1⊙(-(n_Γ⋅∇∇(v)⋅n_Γ) + (γ/h)*∇(v)⋅n_Γ)) * dΓ
               )

        # Define of ghost penalties
        if order >= 2
            # g(u,v) = h^(-2)*(∫( (γg0/h)*jump(u)*jump(v)) * dFg +
            #                   ∫( (γg1*h)*jump(n_Fg⋅∇(u))*jump(n_Fg⋅∇(v)) ) * dFg +
            #                   ∫( (γg2*h^3)*jump_nn(u,n_Fg)*jump_nn(v,n_Fg) ) * dFg)
            g(u,v) = h^(-2)*( ∫( (γg1*h)*jump(n_Fg⋅∇(u))*jump(n_Fg⋅∇(v)) ) * dFg +
                              ∫( (γg2*h^3)*jump_nn(u,n_Fg)*jump_nn(v,n_Fg) ) * dFg)
        else
            println("Not supported order:", order)
        end

        A(u,v) = a(u,v) + g(u,v)


        # Assemble of system
        op = AffineFEOperator(a, l, U, V)
        uh = solve(op)
        A =  get_matrix(op)
        ndof = size(A)[1]
        cond_number = ( 1/sqrt(ndof) )*cond(A,Inf)

        u_inter = interpolate(u_ex, V) # remove?
        e = u_inter - uh
        el2 = sqrt(sum( ∫(e*e)dΩ ))

        # TODO: Add α into ∫(e⊙e)*dΩ
        eh_energy = sqrt(sum( ∫(e⊙e)dΩ + ∫( ∇∇(e)⊙∇∇(e) )dΩ
                      + ( γ/h ) * ∫(jump(∇(e)⋅n_Λ) ⊙ jump(∇(e)⋅n_Λ))dΛ
                      + ( h/γ ) * ∫(mean_nn(e,n_Λ) ⊙ mean_nn(e,n_Λ))dΛ
                      + ( γ/h ) * ∫((∇(e)⋅n_Γ) ⊙ (∇(e)⋅n_Γ))dΓ
                      + ( h/γ ) * ∫(( n_Γ ⋅ ∇∇(e)⋅ n_Γ ) ⊙ ( n_Γ ⋅ ∇∇(e)⋅ n_Γ ))dΓ
                     ))

        eh1 = sqrt(sum( ∫( e⊙e + ∇(e)⊙∇(e) )dΩ ))

        vtkdirname =dirname*"/order_"*string( order )*"_n_"*string(n)
        mkpath(vtkdirname)
        # dirname = vtkdirname * "order_" * string(order) * "_n_" * string(n)


        # Write out models and computational domains for inspection
        writevtk(bgmodel,   vtkdirname*"/model")
        writevtk(Ω,         vtkdirname*"/Omega")
        writevtk(Ω_act,     vtkdirname*"/Omega_act")
        writevtk(Λ,         vtkdirname*"/Lambda")
        writevtk(Γ,         vtkdirname*"/Gamma")
        writevtk(Fg,        vtkdirname*"/Fg")
        writevtk(Λ,         vtkdirname*"/jumps",      cellfields=["jump_u"=>jump(uh)])
        writevtk(Ω,         vtkdirname*"/sol",        cellfields=["e"=>e, "uh"=>uh, "u"=>u_inter])

        sol = Solution( el2=el2, eh1=eh1, eh_energy=eh_energy,
                        cond_number=cond_number, ndof=ndof)
        return sol
    end

end # module


function generate_figures(;ns, el2s, eh1s, ehs_energy, cond_numbers, ndofs, dirname::String)
    filename = dirname*"/conv"

    hs = 1 .// ns
    # hs_str =  latexify.(hs)

    compute_eoc(hs, errs) = log.(errs[1:end-1]./errs[2:end])./log.(hs[1:end-1]./hs[2:end])
    eoc_l2 = compute_eoc(hs, el2s)
    eoc_eh1 = compute_eoc(hs,eh1s)
    eoc_eh_energy = compute_eoc(hs,ehs_energy)

    eoc_l2 =  [nothing; eoc_l2]
    eoc_eh1 =  [nothing; eoc_eh1]
    eoc_eh_energy =  [nothing; eoc_eh_energy]
    header = [L"$n$", L"$\Vert e \Vert_{L^2}$", "EOC", L"$ \Vert e \Vert_{H^1}$", "EOC", L"$\Vert e \Vert_{ a_h,* }$", "EOC", "Cond number", "ndofs"]
    data = hcat(ns, el2s,  eoc_l2, eh1s, eoc_eh1, ehs_energy, eoc_eh_energy, cond_numbers, ndofs)

    formatters = (ft_printf("%.0f", [1]), ft_printf("%.2f", [3, 5, 7]), ft_printf("%.1E", [2, 4, 6, 8, 9]), ft_nonothing)

    open(filename*".tex", "w") do io
        pretty_table(io, data, header=header, backend=Val(:latex ), formatters = formatters )
    end

    minimal_header = ["n", "L2", "EOC", "H1", "EOC", "a_h", "EOC", "cond", "const", "ndofs"]
    data = hcat(ns, el2s,  eoc_l2, eh1s, eoc_eh1, ehs_energy, eoc_eh_energy, cond_numbers, cond_numbers.*hs.^4, ndofs)
    formatters = ( ft_printf("%.0f",[1,10]), ft_printf("%.2f",[3,5,7]), ft_printf("%.1E",[2,4,6,8,9]), ft_nonothing )
    pretty_table(data, header=minimal_header, formatters =formatters )

    # Initial plot with the first data series
    p = Plots.plot(hs, el2s, label=L"\Vert e \Vert_{L^2}", legend=:bottomright, xscale=:log2, yscale=:log2, minorgrid=true)
    Plots.scatter!(p, hs, el2s, primary=false)

    # Add the second data series
    Plots.plot!(p, hs, eh1s, label=L"\Vert e \Vert_{H^1}")
    Plots.scatter!(p, hs, eh1s, primary=false)

    # Add the third data series
    Plots.plot!(p, hs, ehs_energy, label=L"\Vert e \Vert_{a_{h,*}}")
    Plots.scatter!(p, hs, ehs_energy, primary=false)

    # Configs
    Plots.xlabel!(p, "h")
    Plots.plot!(p, xscale=:log2, yscale=:log2, minorgrid=true)
    Plots.plot!(p, legendfontsize=14)  # Adjust the value 12 to your desired font size


    # Save the plot as a .png file using the GR backend
    # Plots.gr()
    Plots.pgfplotsx()
    Plots.savefig(p, filename*"_plot.png")
    Plots.savefig(p, filename*"_plot.tex")
end

function convergence_analysis(;order,  ns, dirname, u_ex)

    el2s = Float64[]
    eh1s = Float64[]
    ehs_energy = Float64[]
    cond_numbers = Float64[]
    ndofs = Float64[]

    for n in ns

        sol = Solver.run( order=order, n=n, u_ex=u_ex, dirname=dirname)

        push!(el2s, sol.el2)
        push!(eh1s, sol.eh1)
        push!(ehs_energy, sol.eh_energy)
        push!(cond_numbers, sol.cond_number)
        push!(ndofs, sol.ndof)
    end
    generate_figures(ns=ns, el2s=el2s, eh1s=eh1s, ehs_energy=ehs_energy,
                     cond_numbers=cond_numbers, ndofs=ndofs, dirname=dirname)
end


function main()

    # %% Manufactured solution
    L, m, r = (1, 1, 1)

    # Examples that works badly
    # u_ex(x) = sin(2π*x[1])*cos(2π*x[2])
    # u_ex(x) = 100*sin(m*( 2π/L )*x[1])*cos(r*( 2π/L )*x[2])
    # u_ex(x) = 100*cos(m*( 2π/L )*x[1])*cos(r*( 2π/L )*x[2])

    # Good examples
    u_ex(x) = (x[1]^2 + x[2]^2  - 1)^2*sin(2π*x[1])*cos(2π*x[2])
    # u_ex(x) = (x[1]^2 + x[2]^2  - 1)^2

    resultdir= "figures/biharmonic_CutCIP/"*string(Dates.now())
    println(resultdir)
    mkpath(resultdir)

    ns = [2^3, 2^4, 2^5, 2^6, 2^7, 2^8]

    @time convergence_analysis( order=2, ns=ns,  dirname=resultdir, u_ex=u_ex,)
end


main()
