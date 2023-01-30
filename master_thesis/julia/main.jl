include("results.jl")
include("biharmonic_equation.jl")


function makedir(dirname)
    if (isdir(dirname))
        rm(dirname, recursive=true)
    end
    mkdir(dirname)
end



function convergence_analysis(; L, m, r, orders, ns, dirname)
    println("Run convergence",)

    for i in 1:length(orders)
        order = orders[i]

        el2s = Float64[]
        eh1s = Float64[]
        ehs_energy = Float64[]
        println("Run convergence tests: order = "*string(order))

        for n in ns
            settings = BiharmonicEquation.SolverSettings(order=order, L=L, n=n, m=m, r=r)
            res = BiharmonicEquation.run_CP_method(settings)
            push!(el2s, res.el2)
            push!(eh1s, res.eh1)
            push!(ehs_energy, res.eh_energy)
        end
        Results.generate_figures(ns=ns, el2s=el2s, eh1s=eh1s, ehs_energy=ehs_energy, order=order, dirname=dirname)
    end

end

function main()
    resultdir= "figures/biharmonic/"*string(Dates.now())
    mkpath(resultdir)

    function run(;  L,m,r)
        orders=[2,3,4]
        ns = [2^2, 2^3, 2^4]#, 2^5]#, 2^6, 2^7]
        dirname = resultdir*"/L_"*string(round(L,digits=2))*"_m_"*string(m)*"_r_"*string(r);
        makedir(dirname)
        convergence_analysis( L=L, m=m, r=r, orders=orders, ns=ns, dirname=dirname)
    end

    run(L=1,m=1,r=1)
end

@time main()
