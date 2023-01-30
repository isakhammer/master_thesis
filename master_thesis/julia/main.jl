include("biharmonic_equation.jl")
using Test
using Plots
import Dates
import CairoMakie
using LaTeXStrings
using Latexify
using PrettyTables

function generate_figures(hs, hs_str, el2s, eh1s, ehs_energy, γ, order::Integer, dirname::String)
    filename = dirname*"/conv_order_"*string(order)*"_gamma_"*string(round( γ,digits=2 ))

    function generate_plot(hs, el2s, eh1s, ehs_energy)
        fig = CairoMakie.Figure()
        ax = CairoMakie.Axis(fig[1, 1], yscale = log10, xscale= log2,
                             yminorticksvisible = true, yminorgridvisible = true, yminorticks = CairoMakie.IntervalsBetween(8),
                             xlabel = L"h/{L}", ylabel = "error norms")

        CairoMakie.lines!(hs, el2s, label= L"$L^2$ norm", linewidth=2)
        CairoMakie.lines!(hs, eh1s, label= L"$H_1$ norm", linewidth=2)
        CairoMakie.lines!(hs, ehs_energy, label= L"$energy$  norm ", linewidth=2)
        CairoMakie.scatter!(hs, el2s)
        CairoMakie.scatter!(hs, eh1s)
        CairoMakie.scatter!(hs, ehs_energy)
        file = filename*".png"
        CairoMakie.Legend(fig[1,2], ax, framevisible = true)
        CairoMakie.save(file,fig)
    end

    function compute_eoc(hs, errs)
        eoc = log.(errs[1:end-1]./errs[2:end])./log.(hs[1:end-1]./hs[2:end])
        return eoc
    end

    function generate_table(;hs_str, hs, el2s, eh1s, ehs_energy)
        eoc_l2 = compute_eoc(hs, el2s)
        eoc_eh1 = compute_eoc(hs,eh1s)
        eoc_eh_energy = compute_eoc(hs,ehs_energy)
        eoc_l2 =  [nothing; eoc_l2]
        eoc_eh1 =  [nothing; eoc_eh1]
        eoc_eh_energy =  [nothing; eoc_eh_energy]

        data = hcat(hs_str, el2s,  eoc_l2, eh1s, eoc_eh1, ehs_energy, eoc_eh_energy)
        header = [L"h/{L} ", L"$L^2$ norm", "EOC", L"$H_1$ norm", "EOC", "energy norm", "EOC"]

        open(filename*".tex", "w") do io
            pretty_table(io, data, header=header, backend=Val(:latex ), formatters = ( ft_printf("%.3E"), ft_nonothing ))
        end

        open(filename*".txt", "w") do io
            pretty_table(io, data, header=header, formatters = ( ft_printf("%.3E"), ft_nonothing ))
        end
    end

    generate_plot(hs, el2s, eh1s, ehs_energy)
    generate_table(hs_str=hs_str, hs=hs, el2s=el2s, eh1s=eh1s, ehs_energy=ehs_energy)
end


function makedir(dirname)
    if (isdir(dirname))
        rm(dirname, recursive=true)
    end
    mkdir(dirname)
end



function convergence_analysis(;figdir, L, u::Function, orders, γs , ns)
    println("Run convergence",)

    hs = 1 .// ns # does not render nice in latex table if L=2π
    hs_str =  latexify.(hs)

    @test length( orders ) == length(γs)
    @test length( ns ) == length(hs)

    for i in 1:length(orders)
        order = orders[i]
        γ = γs[i]

        el2s = Float64[]
        eh1s = Float64[]
        ehs_energy = Float64[]
        println("Run convergence tests: order = "*string(order))

        for n in ns
            res = BiharmonicEquation.run_CP_method(n=n, L=L, γ=γ, order=order, u=u)
            push!(el2s, res.el2)
            push!(eh1s, res.eh1)
            push!(ehs_energy, res.eh_energy)
        end
        generate_figures(hs, hs_str, el2s, eh1s, ehs_energy, γ, order, figdir)
    end

end

# function run_gamma_analysis(;figdir, L, u::Function, orders, γs, ns)
#     hs = 1 .// ns
#     for order in orders
#         println("Run gamma = ", order, " of ", orders)
#         fig = CairoMakie.Figure()
#         ax = CairoMakie.Axis(fig[1, 1], yscale = log10, xscale= log2,
#                              yminorticksvisible = true, yminorgridvisible = true, yminorticks = CairoMakie.IntervalsBetween(8),
#                              xlabel = L"h/{L}", ylabel = L"$L^2$ error" , title="Order: "*string(order))
#         for i in 1:length(γs)
#             γ = γs[i]
#             el2s = Float64[]
#             println("γ ", ": ", i, "/", length(γs))

#             for n in ns
#                 res = BiharmonicEquation.run_CP_method(n=n, L=L, γ=γ, order=order, u=u)
#                 push!(el2s, res.el2)
#             end

#             CairoMakie.lines!(hs, el2s, label=string(γ), linewidth=2)
#             CairoMakie.scatter!(hs, el2s)
#         end

#         CairoMakie.Legend(fig[1,2], ax,L"$\gamma$ values", framevisible = true)
#         CairoMakie.save(figdir*"/gamma_analysis_order"*string(order)*".png",fig)
#     end
# end


function main()
    mainfigdir = "figures"
    if !(isdir(mainfigdir))
        mkdir(mainfigdir)
    end

    # makedir(mainfigdir)
    mainfigdir= mainfigdir*"/"*string(Dates.now())
    makedir(mainfigdir)

    function run(;  L,m,r)
        orders=[2,3,4]
        # ns = [2^2,2^3,2^4,2^5,2^6,2^7,2^8]
        γs_orders  = @. 1.5*orders*( orders+1)

        ns = [2^2, 2^3, 2^4, 2^5]#, 2^6, 2^7]
        γs = [2^0, 2^1, 2^2, 2^3, 2^4, 2^5,2^6]
        figdir = mainfigdir*"/L_"*string(round(L,digits=2))*"_m_"*string(m)*"_r_"*string(r);
        makedir(figdir)
        u = BiharmonicEquation.man_sol(L=L,m=m,r=r)
        convergence_analysis(figdir=figdir, L=L, u=u,  orders=orders, γs=γs_orders, ns=ns)
        # run_gamma_analysis(figdir=figdir, L=L,u=u, orders=orders, γs=γs, ns=ns)
    end

    run(L=1,m=1,r=1)
end

@time main()
