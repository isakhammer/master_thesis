include("biharmonic_equation.jl")
using Test
using Plots
using LaTeXStrings
using PrettyTables

function generate_figures(hs, hs_str, el2s, eh1s, γ::Integer, order::Integer, dirname::String)
    filename = dirname*"/convergence_d_"*string(order)*"_gamma_"*string(γ)

    function slope(hs,errors)
      x = log10.(hs)
      y = log10.(errors)
      linreg = hcat(fill!(similar(x), 1), x) \ y
      linreg[2]
    end

    function generate_plot(hs,el2s, eh1s )
        p_L2 = slope(hs,el2s)
        p_H1 = slope(hs,eh1s)
        println("p_H1: ", p_H1, ", p_L2: ", p_L2)

        p = Plots.plot(hs,[el2s eh1s ],
            xaxis=:log, yaxis=:log,
            label=[L"Error norm in $L_2(\Omega)$ where $p_1 = $"*string(round(p_L2,digits=2)) L"Error norm in $ H^1(\Omega)$ where $p_2 =$"*string(round(p_H1,digits=2)) ],
            # label=[L"Error norm in $L_2(\Omega)$ where $p_1 = $" "1"],
            shape=:auto,
            legend=:topleft,
            legendfontsize=10,
            xlabel=L"$h$",ylabel="error norm" , show = true)

        file = filename*".png"
        Plots.savefig(p,file )
    end

    function generate_table(hs_str::Vector, eh1::Vector{Float64}, el2::Vector{Float64})
        lgl2 = log.(el2[2:end]./el2[1:end-1])
        lgh1 = log.(eh1[2:end]./eh1[1:end-1])
        lgl2 =  [nothing; lgl2]
        lgh1 =  [nothing; lgh1]

        data = hcat(hs_str, el2, eh1, lgl2, lgh1)
        header = ["h", L"$L_2$", L"$H^1$", L"$log_2(e^{2h}_{L^2(\Omega )}/e^{h}_{L^2(\Omega )}) $", L"$log_2(e^{2h}_{H_1(\Omega )}/e^{h}_{H_1(\Omega )}) $"]


        pretty_table(data, header=header, formatters = ( ft_printf("%.3E"), ft_nonothing )) # remove
        # open(filename*".tex", "w") do io
        #     pretty_table(io, data, header=header, backend=Val(:latex ), formatters = ( ft_printf("%.3E"), ft_nonothing ))
        # end
    end

    println("Generating plot", hs ,el2s)
    generate_plot(hs, el2s, eh1s)

    println("Generating table")
    generate_table(hs_str, el2s, eh1s)
end


# Generate plots
function makedir(dirname)
    if (isdir(dirname))
        rm(dirname, recursive=true)
    end
    mkdir(dirname)
end

function run_examples(dirname)
    makedir(dirname)
    ndir(n) = dirname*"/n_"*string(n)
    ns = [100]
    order=2
    for n in ns
        ss = BiharmonicEquation.generate_square_spaces(n=n, L=2π, order=order)
        sol = BiharmonicEquation.run_CP_method(ss=ss)
        BiharmonicEquation.generate_vtk(ss=ss,sol=sol,dirname=ndir(n))
        @test sol.el2 < 10^-1
    end
end



function convergence_analysis(dirname)
    orders = [2,3,4]
    γs = [1, 20, 40]
    ns  = [2^3,2^4,2^5,2^6,2^7]
    hs = 1 .// ns # does not work for 2π
    hs_str =  latexify.(hs)

    for i in 1:length(orders)
        order = orders[i]
        γ = γs[i]

        el2s = Float64[]
        eh1s = Float64[]
        println("Run convergence tests: order = "*string(order))

        for n in ns
            ss = BiharmonicEquation.generate_square_spaces(n=n,L=2π, γ=γ, order=order)
            sol = BiharmonicEquation.run_CP_method(ss=ss)
            push!(el2s, sol.el2)
            push!(eh1s, sol.eh1)
            # push!(hs,   ss.h)
        end

        println("Generate figures")
        generate_figures(hs, hs_str, el2s, eh1s, γ, order, dirname)
    end

end

function main()

    folder = "biharmonic_julia_test_results"
    exampledir = folder*"/example"
    makedir(folder)
    run_examples(exampledir)

    println("Generating figures")
    figdir = "figures"
    makedir(figdir)
    convergence_analysis(figdir)
end

main()
