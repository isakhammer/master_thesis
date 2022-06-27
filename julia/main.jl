include("biharmonic_equation.jl")
using Test
using Plots
import Dates
import CairoMakie
using LaTeXStrings
using Latexify
using PrettyTables

function generate_figures(hs, hs_str, el2s, ehs, γ::Integer, order::Integer, dirname::String)
    filename = dirname*"/conv_order_"*string(order)*"_gamma_"*string(γ)

    function generate_plot(hs,el2s, ehs )
        fig = CairoMakie.Figure()
        ax = CairoMakie.Axis(fig[1, 1], yscale = log10, xscale= log2,
                             yminorticksvisible = true, yminorgridvisible = true, yminorticks = CairoMakie.IntervalsBetween(8),
                             xlabel = "h", ylabel = "error norms")

        CairoMakie.lines!(hs, el2s, label= L"L2 norm", linewidth=2)
        CairoMakie.lines!(hs, ehs, label= L"energy norm", linewidth=2)
        CairoMakie.scatter!(hs, el2s)
        CairoMakie.scatter!(hs, ehs)
        file = filename*".png"
        CairoMakie.Legend(fig[1,2], ax, framevisible = true)
        CairoMakie.save(file,fig)
    end

    function compute_eoc(hs, errs)
        eoc = log.(errs[1:end-1]./errs[2:end])./log.(hs[1:end-1]./hs[2:end])
        return eoc
    end

    function generate_table(;hs_str, hs, el2s, ehs)
        eoc_l2 = compute_eoc(hs, el2s)
        eoc_eh = compute_eoc(hs,ehs)
        eoc_l2 =  [nothing; eoc_l2]
        eoc_eh =  [nothing; eoc_eh]

        data = hcat(hs_str, el2s,  eoc_l2, ehs, eoc_eh)
        header = [L"h", L"$L_2$ norm", "EOC", "energy norm", "EOC"]

        open(filename*".tex", "w") do io
            pretty_table(io, data, header=header, backend=Val(:latex ), formatters = ( ft_printf("%.3E"), ft_nonothing ))
        end

        open(filename*".txt", "w") do io
            pretty_table(io, data, header=header, formatters = ( ft_printf("%.3E"), ft_nonothing ))
        end
    end

    generate_plot(hs, el2s, ehs)
    generate_table(hs_str=hs_str, hs=hs, el2s=el2s, ehs=ehs)
end


function makedir(dirname)
    if (isdir(dirname))
        rm(dirname, recursive=true)
    end
    mkdir(dirname)
end

function run_examples(;figdir, L, u::Function, ns = [2^3,2^5], γ=2, order=2)
    println("Run Examples")
    exampledir = figdir*"/example"
    makedir(exampledir)
    ndir(n) = exampledir*"/n_"*string(n)

    for n in ns
        res = BiharmonicEquation.run_CP_method(n=n, L=L, γ=γ, order=order)
        BiharmonicEquation.generate_vtk(res=res,dirname=ndir(n))
        # @test sol.el2 < 10^0
    end
end


function convergence_analysis(;figdir, L, u::Function, orders = [2,3,4], γs = [2, 8, 16], ns = [2^3,2^4,2^5,2^6,2^7])
    println("Run convergence",)

    hs = 1 .// ns # does render nice in latex table if L=2π
    hs_str =  latexify.(hs)

    @test length( orders ) == length(γs)
    @test length( ns ) == length(hs)

    for i in 1:length(orders)
        order = orders[i]
        γ = γs[i]

        el2s = Float64[]
        ehs = Float64[]
        println("Run convergence tests: order = "*string(order))

        for n in ns
            res = BiharmonicEquation.run_CP_method(n=n, L=L, γ=γ, order=order, u=u)
            push!(el2s, res.el2)
            push!(ehs, res.eh)
            # push!(hs,   ss.h)
        end
        generate_figures(hs, hs_str, el2s, ehs, γ, order, figdir)
    end

end

function run_gamma_analysis(;figdir, L, u::Function, orders = [2,3,4], γs = [2^0,2^1, 2^2,2^3,2^4,2^5,2^6, 2^7], ns = [2^3,2^4,2^5,2^6,2^7])
    hs = 1 .// ns
    for order in orders
        println("Run gamma = ", order, " of ", orders)
        fig = CairoMakie.Figure()
        ax = CairoMakie.Axis(fig[1, 1], yscale = log10, xscale= log2,
                             yminorticksvisible = true, yminorgridvisible = true, yminorticks = CairoMakie.IntervalsBetween(8),
                             xlabel = "h", ylabel = "L2-error" , title="Order: "*string(order))
        for i in 1:length(γs)
            γ = γs[i]
            el2s = Float64[]
            ehs = Float64[]
            println("γ ", ": ", i, "/", length(γs))

            for n in ns
                res = BiharmonicEquation.run_CP_method(n=n, L=L, γ=γ, order=order, u=u)
                push!(el2s, res.el2)
                push!(ehs, res.eh)
                # push!(hs,   ss.h)
            end
            CairoMakie.lines!(hs, el2s, label=string(γ), linewidth=2)
            CairoMakie.scatter!(hs, el2s)
        end

        CairoMakie.Legend(fig[1,2], ax,L"$\gamma$ values", framevisible = true)
        CairoMakie.save(figdir*"/gamma_analysis_order"*string(order)*".png",fig)
    end
end


function main()
    mainfigdir = "figures"

    if !(isdir(mainfigdir))
        mkdir(mainfigdir)
    end

    # makedir(mainfigdir)
    mainfigdir= mainfigdir*"/"*string(Dates.now())
    makedir(mainfigdir)

    function run(;  L,m,r, orders=[2,3,4], γs=[2,8,32])
        figdir = mainfigdir*"/L_"*string(round(L,digits=2))*"_m_"*string(m)*"_r_"*string(r);
        makedir(figdir)
        u = BiharmonicEquation.man_sol(L=L,m=m,r=r)
        convergence_analysis(figdir=figdir, L=L, u=u,  orders=orders, γs=γs)
        run_gamma_analysis(figdir=figdir, L=L,u=u)
        # run_examples(figdir=figdir, L=L,u=u)
    end

    run(L=2π, m=1,r=1, orders=[2,3,4], γs=[2,8,16])
    run(L=1,m=1,r=1, orders=[2,3,4], γs=[2,8,16])
    run(L=1,m=7,r=3, orders=[2,3,4], γs=[2,8,16])
end


@time main()
