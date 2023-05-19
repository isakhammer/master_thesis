include("cahn_hilliard_CutCIP.jl")
include("nonlinear_cahn_hilliard_CutCIP.jl")
include("heat_eq.jl")

using Test
using Plots
using PrettyTables
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

    if (false)
        # # Create a matrix to store the error values for each combination
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
    end

    # Function to compute a corresponding EOC matrix
    hs = 1 .// ns


    compute_eoc_vector(hs, dts, error_vector) = log.(error_vector[1:end-1]./error_vector[2:end])./( log.(hs[1:end-1]./hs[2:end]) + log.(dts[1:end-1]./dts[2:end]) ) # Error! It only works for vectors

    function compute_eoc_matrices(error_matrix, hs, dts)
        N,M = size(error_matrix)

        transient_eoc_matrix = Matrix{Union{Nothing, Float64}}(undef, N, M)
        # Transient EOC
        for i in 1:N
            error_vector = error_matrix[i,:]
            his = fill(hs[i], N)
            eoc_vector = compute_eoc_vector(his, dts, error_vector)
            eoc_vector =  [nothing; eoc_vector]
            transient_eoc_matrix[i,:] = eoc_vector
        end

        spatial_eoc_matrix = Matrix{Union{Nothing, Float64}}(undef, N, M)
        # Transient EOC
        for j in 1:M
            error_vector = error_matrix[:,j]
            dtjs = fill(dts[j], M)
            eoc_vector = compute_eoc_vector(hs, dtjs, error_vector)
            eoc_vector =  [nothing; eoc_vector]
            spatial_eoc_matrix[:,j] = eoc_vector
        end

        # diagonal EOC
        # for j in 1:M
        #     for i in 1:N
        #         error_vector = error_matrix[:,j]
        #         dtjs = fill(dts[j], M)
        #         eoc_vector = compute_eoc_vector(hs, dtjs, error_vector)
        #         eoc_vector =  [nothing; eoc_vector]
        #         spatial_eoc_matrix[:,j] = eoc_vector
        #     end
        # end

        return transient_eoc_matrix, spatial_eoc_matrix
    end

    # Test example
    error_matrix = ones(length(hs), length(dts))
    pretty_table(error_matrix)
    transient_eoc_matrix, spatial_eoc_matrix = compute_eoc_matrices(error_matrix, hs, dts)

    pretty_table(transient_eoc_matrix)
    pretty_table(spatial_eoc_matrix)


    # Create LaTeX-formatted strings for dts and ns
    dt_str = latexify.(1 .// Int.(1 ./ dts))
    hs_str = latexify.(1 .// ns)

    # Create the header for the table
    header = [" $(hs_str[i])" for i in 1:length(hs_str)]

    # Create a table with row names as the dt values
    row_names = dt_str

    # Print the EOC matrix as a pretty table
    pretty_table(transient_eoc_matrix, header, row_names=row_names, formatters = ft_printf("%12.5f"))
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
