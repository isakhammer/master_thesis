include("biharmonic_AgFEM_CIP_laplace.jl")
# include("biharmonic_CutCIP.jl")

using Dates
import Plots

# default_size = (800, 600)
default_size = (400, 300)

# Plots.pgfplotsx()
# endfix=".tex"
Plots.gr()
endfix=".png"

using LaTeXStrings
using Latexify
using PrettyTables


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

    open(filename*".txt", "w") do io
        pretty_table(io, data, header=minimal_header, formatters=formatters)
    end
    # Initial plot with the first data series
    p = Plots.plot(hs, el2s, label=L"\Vert e \Vert_{L^2}", size=default_size, legend=:outertopright, xscale=:log2, yscale=:log2, minorgrid=true)
    Plots.scatter!(p, hs, el2s, primary=false)

    # Add the second data series
    Plots.plot!(p, hs, eh1s, label=L"\Vert e \Vert_{H^1}")
    Plots.scatter!(p, hs, eh1s, primary=false)

    # Add the third data series
    Plots.plot!(p, hs, ehs_energy, label=L"\Vert e \Vert_{a_{h,*}}")
    Plots.scatter!(p, hs, ehs_energy, primary=false)

    # Configs
    Plots.xlabel!(p, "h")
    Plots.ylabel!(p, L"\Vert e \Vert_{}")
    Plots.plot!(p, xscale=:log2, yscale=:log2, minorgrid=true)
    Plots.plot!(p, legendfontsize=14)  # Adjust the value 12 to your desired font size

    Plots.savefig(p, filename*"_plot"*endfix)
end

function convergence_analysis(; type::SolverLaplace.CutFEMType, order, ns, dirname, u_ex, run_solver::Function, L=1.11, δ=0.0, γ=10.0, γg1=5, γg2=0.01)

    el2s = Float64[]
    eh1s = Float64[]
    ehs_energy = Float64[]
    cond_numbers = Float64[]
    ndofs = Float64[]

    println("Convergence test", ns)
    for n in ns

        sol = run_solver(; type, order, n, u_ex, dirname, L=L, δ=δ, γ=γ, γg1=γg1, γg2=γg2)

        push!(el2s, sol.el2)
        push!(eh1s, sol.eh1)
        push!(ehs_energy, sol.eh_energy)
        push!(cond_numbers, sol.cond_number)
        push!(ndofs, sol.ndof)
    end
    generate_figures(ns=ns, el2s=el2s, eh1s=eh1s, ehs_energy=ehs_energy,
                     cond_numbers=cond_numbers, ndofs=ndofs, dirname=dirname)
end


function main(; type::SolverLaplace.CutFEMType, order, ns, γ, γg1, γg2)
    # %% Manufactured solution
    L, m, r = (1, 1, 1)
    u_ex(x) = (x[1]^2 + x[2]^2 - 1)^2*sin(m*( 2π/L )*x[1])*cos(r*( 2π/L )*x[2])

    # Laplace
    laplace =true
    # hessian =false

    if laplace
        resultdir= "figures/eoc_test/laplace_"*string(Dates.now())
        println(resultdir)
        mkpath(resultdir)
        @time convergence_analysis(type=type, order=order, ns=ns, dirname=resultdir, u_ex=u_ex, run_solver=SolverLaplace.run,  
                                   L=1.12, δ=0.0, γ=γ, γg1=γg1, γg2=γg2)

    end

    # if hessian
    #     resultdir= "figures/eoc_test/hessian_"*string(Dates.now())
    #     println(resultdir)
    #     mkpath(resultdir)
    #     ns = [2^3, 2^4, 2^5, 2^6, 2^7, 2^8]
    #     @time convergence_analysis( ns=ns,  dirname=resultdir, u_ex=u_ex, run_solver=SolverHessian.run,  L=2.11, δ=0.0, γ=10, γg1=5, γg2=0.01)

    # end
end


ns = [2^3, 2^4, 2^5, 2^6, 2^7, 2^8]
# ns = [2^4, 2^5, 2^6, 2^7]
# type=SolverLaplace.AgFEM
type=SolverLaplace.CutFEM
order = 2
γ = 5*order*(order-1)
γg1=γ
γg2=0.01
main(type=type, order=order, ns=ns, γ=γ, γg1=γg1, γg2=γg2)
