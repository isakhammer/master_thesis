include("cahn_hilliard_CutCIP.jl")
include("nonlinear_cahn_hilliard_CutCIP.jl")
include("heat_eq.jl")

using Test
using Plots
using Dates
using LaTeXStrings
using Latexify
using PrettyTables

function run_problem(;problem::String, ode_method::String,
        n::Number, dt::Number, u_ex::Function,
        vtkdirname::String)
    if problem == "CH"
        return CH.run(n=n, dt=dt, u_ex=u_ex, ode_method=ode_method, vtkdirname=vtkdirname)
    elseif problem == "NLCH"
        return NLCH.run(n=n, dt=dt, u_ex=u_ex, ode_method=ode_method, vtkdirname=vtkdirname)
    elseif problem == "HE"
        return HE.run(n=n, dt=dt, u_ex=u_ex, ode_method=ode_method, vtkdirname=vtkdirname)
    else
        error("Invalid problem: $problem")
    end
end

function convergence_matrix(; ns::Vector, dts::Vector,
        main_dirname::String, u_ex::Function, problem::String, ode_method::String)

    # Create a matrix to store the error values for each combination
    el2_L2_matrix         = Array{Float64}(undef, length(dts), length(ns))
    eh1_L2_matrix         = Array{Float64}(undef, length(dts), length(ns))
    eh_energy_L2_matrix   = Array{Float64}(undef, length(dts), length(ns))
    el2_inf_matrix        = Array{Float64}(undef, length(dts), length(ns))
    eh1_inf_matrix        = Array{Float64}(undef, length(dts), length(ns))
    eh_energy_inf_matrix  = Array{Float64}(undef, length(dts), length(ns))

    dirname = main_dirname*"/$problem"*"_$ode_method"
    println(dirname)
    mkpath(dirname)
    println("Run convergence for problem $problem and with ODE solver $ode_method")

    # Fill the error matrices
    for i in 1:length(ns)
        for j in 1:length(dts)
            filename = dirname*"_n_$(ns[i])"*"_dt_$(dts[j])"
            sol = run_problem(problem=problem, ode_method=ode_method, n=ns[i],
                              dt=dts[j], u_ex=u_ex, vtkdirname=filename)
            el2_L2_matrix[i,j] = sol.el2s_L2
            el2_L2_matrix[i, j] = sol.el2s_L2
            eh1_L2_matrix[i, j] = sol.eh1s_L2
            eh_energy_L2_matrix[i, j] = sol.ehs_energy_L2
            el2_inf_matrix[i, j] = sol.el2s_inf
            eh1_inf_matrix[i, j] = sol.eh1s_inf
            eh_energy_inf_matrix[i, j] = sol.ehs_energy_inf
        end
    end

    # Function to compute a corresponding EOC matrix
    hs = 1 .// ns
    compute_eoc(hs, dts, error_vector) = log.(errs[1:end-1]./errs[2:end])./( log.(hs[1:end-1]./hs[2:end]) + log.(dts[1:end-1]./dts[2:end]) ) # Error! It only works for vectors

    # Create LaTeX-formatted strings for dts and ns
    dt_str = latexify.(1 .// Int.(1 ./ dts))
    n_str = latexify.(ns)

    # # Create the header for the table
    # header = ["dt \\ n"]
    # append!(header, ["n = $(n_str[i])" for i in 1:length(ns)])

    # # Create a table with row names as the dt values
    # row_names = dt_str

    # # Print the EOC matrix as a pretty table
    # pretty_table(eoc_matrix, header, row_names=row_names, formatters = ft_printf("%12.5f"))
end

function main_convergence()

    main_dirname= "figures/convergence_analysis_"*string(Dates.now())
    println(main_dirname)
    mkpath(main_dirname)

    u_ex(x,t::Real) = sin(t)*(x[1]^2 + x[2]^2 - 1 )^3*sin(x[1])*cos(x[2])
    u_ex(t::Real) = x -> u_ex(x,t)

    dts = [2^-2, 2^-3, 2^-4, 2^-5]
    ns = [2^3, 2^4, 2^5, 2^6]


    @time convergence_matrix( ns=ns, dts=dts, main_dirname=main_dirname,
                             u_ex=u_ex, problem="HE", ode_method="BE")


end

@time main_convergence()
