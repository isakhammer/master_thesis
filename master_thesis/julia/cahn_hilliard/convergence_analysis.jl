include("cahn_hilliard_CutCIP.jl")
using Test
using Plots
using Dates
using LaTeXStrings
using Latexify
using PrettyTables

function generate_figures(ns::Vector, dt,
        el2s_L2::Vector, eh1s_L2::Vector, ehs_energy_L2::Vector,
        el2s_inf::Vector, eh1s_inf::Vector, ehs_energy_inf::Vector,
        dirname::String)

    dts = fill(dt, length(ns))

    # call the original function with ns instead of n
    generate_figures(ns, dts, el2s_L2, eh1s_L2, ehs_energy_L2,
                     el2s_inf, eh1s_inf, ehs_energy_inf, dirname)
end

function generate_figures(n, dts::Vector,
        el2s_L2::Vector, eh1s_L2::Vector, ehs_energy_L2::Vector,
        el2s_inf::Vector, eh1s_inf::Vector, ehs_energy_inf::Vector,
        dirname::String)

    ns = fill(n, length(dts))

    # call the original function with ns instead of n
    generate_figures(ns, dts, el2s_L2, eh1s_L2, ehs_energy_L2,
                     el2s_inf, eh1s_inf, ehs_energy_inf, dirname)
end

function generate_figures(ns::Vector, dts::Vector,
        el2s_L2::Vector, eh1s_L2::Vector, ehs_energy_L2::Vector,
        el2s_inf::Vector, eh1s_inf::Vector, ehs_energy_inf::Vector,
        filename::String)

    hs = 1 .// ns
    compute_eoc(hs, dts, errs) = log.(errs[1:end-1]./errs[2:end])./( log.(hs[1:end-1]./hs[2:end]) + log.(dts[1:end-1]./dts[2:end]) )

    hs_str =  latexify.(hs)
    dt_str =  latexify.( 1 .//Int.(1 ./ dts))

    eoc_l2s_L2 = compute_eoc(hs, dts, el2s_L2)
    eoc_eh1s_L2 = compute_eoc(hs, dts, eh1s_L2)
    eoc_ehs_energy_L2 = compute_eoc(hs,dts, ehs_energy_L2)
    eoc_l2s_inf = compute_eoc(hs, dts, el2s_inf)
    eoc_eh1s_inf = compute_eoc(hs, dts, eh1s_inf)
    eoc_ehs_energy_inf = compute_eoc(hs, dts, ehs_energy_inf)

    eoc_l2s_L2 =  [nothing; eoc_l2s_L2]
    eoc_eh1s_L2 =  [nothing; eoc_eh1s_L2]
    eoc_ehs_energy_L2 =  [nothing; eoc_ehs_energy_L2]
    eoc_l2s_inf =  [nothing; eoc_l2s_inf]
    eoc_eh1s_inf =  [nothing; eoc_eh1s_inf]
    eoc_ehs_energy_inf =  [nothing; eoc_ehs_energy_inf]

    data = hcat(hs_str, dt_str,
                el2s_L2,  eoc_l2s_L2, eh1s_L2, eoc_eh1s_L2, ehs_energy_L2, eoc_ehs_energy_L2,
                el2s_inf,  eoc_l2s_inf, eh1s_inf, eoc_eh1s_inf, ehs_energy_inf, eoc_ehs_energy_inf)

    formatters = (ft_nonothing, ft_nonothing, ft_printf("%.2f", [2, 4, 6, 8, 10, 12]),
                  ft_printf("%.1E", [3, 5, 7, 9, 11, 13]))

    header = ["hs", "dt",
              "L2L2", "EOC", "L2H1", "EOC", "L2ah", "EOC",
              "infL2", "EOC", "infH1", "EOC", "infah", "EOC"]

    pretty_table(data, header=header, formatters =formatters )


    # formatters = ( ft_nonothing, ft_nonothing, ft_printf("%.5f", [3, 5, 7, 9, 11,13]),
    #               ft_printf("%.1E", [2, 4, 6, 8, 10,12]))

    open(filename*".tex", "w") do io
        pretty_table(io, data, header=header, backend=Val(:latex ), formatters = formatters )
    end


end


function convergence_analysis(; ns::Vector, dts::Vector, dirname::String, solver_config, spatial=false, dt_const=2^-3, transient=false, n_const=2^4, diagonal=false)
    println("Run convergence",)

    if (transient)
        el2s_L2 = Float64[]
        eh1s_L2 = Float64[]
        ehs_energy_L2 = Float64[]
        el2s_inf = Float64[]
        eh1s_inf = Float64[]
        ehs_energy_inf = Float64[]
        # Transient EOC
        println("Run transient EOC tests with constant n = "*string(n_const))

        filename = dirname*"/conv_transient"
        for dt in dts
            sol = Solver.run(n=n_const, dt=dt, solver_config=solver_config, vtkdirname=filename)

            push!(el2s_L2, sol.el2s_L2)
            push!(eh1s_L2, sol.eh1s_L2)
            push!(ehs_energy_L2, sol.ehs_energy_L2)
            push!(el2s_inf, sol.el2s_inf)
            push!(eh1s_inf, sol.eh1s_inf)
            push!(ehs_energy_inf, sol.ehs_energy_inf)
        end
        generate_figures(n_const, dts,
                         el2s_L2, eh1s_L2, ehs_energy_L2,
                         el2s_inf, eh1s_inf, ehs_energy_inf,
                         filename)
    end

    if (spatial)
        el2s_L2 = Float64[]
        eh1s_L2 = Float64[]
        ehs_energy_L2 = Float64[]
        el2s_inf = Float64[]
        eh1s_inf = Float64[]
        ehs_energy_inf = Float64[]
        # Spatial EOC
        println("Run spatial EOC tests with constant dt = "*string(dt_const))
        filename = dirname*"/conv_spatial"
        for n in ns
            sol = Solver.run(n=n, dt=dt_const, solver_config=solver_config, vtkdirname=filename)
            push!(el2s_L2, sol.el2s_L2)
            push!(eh1s_L2, sol.eh1s_L2)
            push!(ehs_energy_L2, sol.ehs_energy_L2)
            push!(el2s_inf, sol.el2s_inf)
            push!(eh1s_inf, sol.eh1s_inf)
            push!(ehs_energy_inf, sol.ehs_energy_inf)
        end

        generate_figures(ns, dt_const,
                         el2s_L2, eh1s_L2, ehs_energy_L2,
                         el2s_inf, eh1s_inf, ehs_energy_inf,
                         filename)
    end


    if (diagonal)
        el2s_L2 = Float64[]
        eh1s_L2 = Float64[]
        ehs_energy_L2 = Float64[]
        el2s_inf = Float64[]
        eh1s_inf = Float64[]
        ehs_energy_inf = Float64[]

        # Spatial EOC
        println("Run convergence tests with dt = $dts"*" and ns= $ns")

        if length(ns) != length(dts)
            error("Cannot compute diagonal. Length does not match")
        end

        filename = dirname*"/conv_diagonal"
        for i in 1:length(ns)
            ni = ns[i]
            dti = dts[i]
            sol = Solver.run(n=ni, dt=dti, solver_config=solver_config, vtkdirname=filename)
            push!(el2s_L2, sol.el2s_L2)
            push!(eh1s_L2, sol.eh1s_L2)
            push!(ehs_energy_L2, sol.ehs_energy_L2)
            push!(el2s_inf, sol.el2s_inf)
            push!(eh1s_inf, sol.eh1s_inf)
            push!(ehs_energy_inf, sol.ehs_energy_inf)
        end

        generate_figures(ns, dts,
                         el2s_L2, eh1s_L2, ehs_energy_L2,
                         el2s_inf, eh1s_inf, ehs_energy_inf,
                         filename)
    end

end

function convergence_matrix(; ns::Vector, dts::Vector, dirname::String, solver_config)

    # Create a matrix to store the error values for each combination
    el2_L2_matrix         = Array{Float64}(undef, length(dts), length(ns))
    eh1_L2_matrix         = Array{Float64}(undef, length(dts), length(ns))
    eh_energy_L2_matrix   = Array{Float64}(undef, length(dts), length(ns))
    el2_inf_matrix        = Array{Float64}(undef, length(dts), length(ns))
    eh1_inf_matrix        = Array{Float64}(undef, length(dts), length(ns))
    eh_energy_inf_matrix  = Array{Float64}(undef, length(dts), length(ns))

    # Fill the error matrices
    for i in 1:length(ns)
        for j in 1:length(dts)
            sol = Solver.run(n=n[i], dt=dt[j],
                             solver_config=solver_config,
                             vtkdirname=dirname*"/conv_matrix")
            el2_L2_matrix[i,j] = sol.el2s_L2
            el2_L2_matrix[i, j] = sol.el2s_L2
            eh1_L2_matrix[i, j] = sol.eh1s_L2
            eh_energy_L2_matrix[i, j] = sol.eh_energys_L2
            el2_inf_matrix[i, j] = sol.el2s_inf
            eh1_inf_matrix[i, j] = sol.eh1s_inf
            eh_energy_inf_matrix[i, j] = sol.eh_energys_inf
        end
    end

    # Function to compute a corresponding EOC matrix
    hs = 1 .// ns
    compute_eoc(hs, dts, error_vector) = log.(errs[1:end-1]./errs[2:end])./( log.(hs[1:end-1]./hs[2:end]) + log.(dts[1:end-1]./dts[2:end]) ) # Error! It only works for vectors

    # Create LaTeX-formatted strings for dts and ns
    dt_str = latexify.(1 .// Int.(1 ./ dts))
    n_str = latexify.(ns)

    # Create the header for the table
    header = ["dt \\ n"]
    append!(header, ["n = $(n_str[i])" for i in 1:length(ns)])

    # Create a table with row names as the dt values
    row_names = dt_str

    # Print the EOC matrix as a pretty table
    pretty_table(eoc_matrix, header, row_names=row_names, formatters = ft_printf("%12.5f"))
end

function main_convergence()

    dirname= "figures/cahn_hilliard_CutCIP/example"*string(Dates.now())
    println(dirname)
    mkpath(dirname)

    u_ex(x,t::Real) = sin(t)*(x[1]^2 + x[2]^2 - 1 )^3*sin(x[1])*cos(x[2])
    u_ex(t::Real) = x -> u_ex(x,t)
    exact_sol = Solver.man_sol(u_ex)
    circle = true
    solver_config = Solver.Config(exact_sol, circle)

    dts = [2^-2,2^-3,2^-4,2^-5,2^-6]
    ns = [2^4,2^5,2^6,2^7, 2^8]
    @time convergence_analysis( ns=ns, dts=dts, dirname=dirname, solver_config=solver_config, spatial=true, dt_const=2^-6, transient=true, n_const=2^8, diagonal=true)

end

main_convergence()
