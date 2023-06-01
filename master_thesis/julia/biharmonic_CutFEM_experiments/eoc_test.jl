include("biharmonic_CutCIP_laplace.jl")
include("biharmonic_CutCIP_hessian.jl")
using Statistics

using Dates
using Test
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
using DataFrames
using CSV

# Function to compute EOC and prepend nothing for each pair of hs and errs
function compute_eoc(hs::Vector, errs::Vector)
    eoc = log.(errs[1:end-1] ./ errs[2:end]) ./ log.(hs[1:end-1] ./ hs[2:end])
    return [NaN; eoc]
end

# Function to compute EOC for three pairs of hs and errs vectors
function compute_eoc(hs::Vector, errs1::Vector, errs2::Vector, errs3::Vector)
    eoc1 = compute_eoc(hs, errs1)
    eoc2 = compute_eoc(hs, errs2)
    eoc3 = compute_eoc(hs, errs3)
    return eoc1, eoc2, eoc3
end


function generate_figures(;ns, el2s, eh1s, ehs_energy, cond_numbers, ndofs, dirname::String)
    filename = dirname*"/conv"

    hs = 1 .// ns

    hs_str = ["1/$(n)" for n in ns]

    eoc_l2, eoc_eh1, eoc_eh_energy = compute_eoc(hs, el2s, eh1s, ehs_energy)

    data = hcat(hs_str, ns, el2s,  eoc_l2, eh1s, eoc_eh1, ehs_energy, eoc_eh_energy, cond_numbers, ndofs)
    formatters = (ft_printf("%s", [1]), ft_printf("%.0f", [2]), ft_printf("%.2f", [4, 6, 8]), ft_printf("%.1E", [3, 5, 7, 9, 10]), ft_nonothing)
    minimal_header = ["h","n", "L2", "EOC", "H1", "EOC", "a_h", "EOC", "cond",  "ndofs"]
    pretty_table(data, header=minimal_header, formatters =formatters )
    open(filename*".txt", "w") do io
        pretty_table(io, data, header=minimal_header, formatters=formatters)
    end

    # Producing a CSV file
    CSV.write(filename*".csv", DataFrame(ns=ns, el2s=el2s, eh1s=eh1s, ehs_energy=ehs_energy,
                                         cond_numbers=cond_numbers, ndofs=ndofs), delim=',')

end

function convergence_analysis(; ns, dirname, u_ex, run_solver::Function, L=1.11, δ=0.0, γ=20, γg1=10, γg2=0.01, geometry::String="circle")

    el2s = Float64[]
    eh1s = Float64[]
    ehs_energy = Float64[]
    cond_numbers = Float64[]
    ndofs = Float64[]

    println("Convergence test n = $ns, geo = $geometry" )
    for n in ns

        sol, _ = run_solver(; n, u_ex, dirname, L=L, δ=δ, γ=γ, γg1=γg1, γg2=γg2, geometry=geometry)

        push!(el2s, sol.el2)
        push!(eh1s, sol.eh1)
        push!(ehs_energy, sol.eh_energy)
        push!(cond_numbers, sol.cond_number)
        push!(ndofs, sol.ndof)
    end
    generate_figures(ns=ns, el2s=el2s, eh1s=eh1s, ehs_energy=ehs_energy,
                     cond_numbers=cond_numbers, ndofs=ndofs, dirname=dirname)

    hs = 1 .//ns

    @testset "EOC tests" begin
        expected_eoc_l2, expected_eoc_h1, expected_eoc_energy = (2, 2, 1)

        # Compute the mean of all EOC except for the "nothing" value and first one.
        eoc_l2, eoc_eh1, eoc_eh_energy = compute_eoc(hs, el2s, eh1s, ehs_energy)
        mean_eoc_l2 = mean(eoc_l2[3:end])
        mean_eoc_h1 = mean(eoc_eh1[3:end])
        mean_eoc_energy = mean(eoc_eh_energy[3:end])

        # Check if the computed means are approximately equal to the expected values
        @test isapprox(mean_eoc_l2, expected_eoc_l2, atol=0.25)
        @test isapprox(mean_eoc_h1, expected_eoc_h1, atol=0.25)
        @test isapprox(mean_eoc_energy, expected_eoc_energy, atol=0.25)
    end

end


function main()


    maindir = "figures/eoc-test"
    if isdir(maindir)
        rm(maindir; recursive=true)
        mkdir(maindir)
    end
    # Manufactured solution
    l, m, r = (2, 1, 1)
    u_ex(x) = sin(m*( 2π/l )*x[1])*cos(r*( 2π/l )*x[2])

    γ, γg1, γg2 = 20, 10, 1
    L, δ = 2.7, 0.0
    ns = [2^3, 2^4, 2^5, 2^6, 2^7, 2^8, 2^9]

    @testset "Laplace Flower EOC tests" begin
        geometry ="flower"
        resultdir= "$maindir/eoc-laplace-$(geometry)-L-$(L)-gamma0-$(γ)-gamma1-$(γg1)-gamma2-$(γg2)"
        println(resultdir)
        mkpath(resultdir)
        @time convergence_analysis( ns=ns,  dirname=resultdir, u_ex=u_ex, run_solver=SolverLaplace.run,  L=L, δ=δ, γ=γ, γg1=γg1, γg2=γg2, geometry=geometry)
    end

    # Manufactured solution
    u_ex(x) = (x[1]^2 + x[2]^2 - 1)^2*sin(m*( 2π/l )*x[1])*cos(r*( 2π/l )*x[2])
    @testset "Laplace EOC tests" begin
        geometry ="circle"
        resultdir= "$maindir/eoc-laplace-$(geometry)-L-$(L)-gamma0-$(γ)-gamma1-$(γg1)-gamma2-$(γg2)"
        println(resultdir)
        mkpath(resultdir)
        @time convergence_analysis( ns=ns,  dirname=resultdir, u_ex=u_ex, run_solver=SolverLaplace.run,  L=L, δ=δ, γ=γ, γg1=γg1, γg2=γg2, geometry=geometry)
    end

    @testset "Hessian EOC tests" begin
        geometry ="circle"
        resultdir= "$maindir/eoc-hessian-$(geometry)-L-$(L)-gamma0-$(γ)-gamma1-$(γg1)-gamma2-$(γg2)"
        println(resultdir)
        mkpath(resultdir)
        @time convergence_analysis( ns=ns,  dirname=resultdir, u_ex=u_ex, run_solver=SolverHessian.run,  L=L, δ=δ, γ=γ, γg1=γg1, γg2=γg2)
    end
end

@testset "Main" begin
    main()
end
