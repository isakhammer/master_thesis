
using PrettyTables
using LaTeXStrings

function generate_table(h::Vector{String}, eh1::Vector{Float64}, el2::Vector{Float64})
    lgl2 = log.(el2[2:end]./el2[1:end-1])
    lgh1 = log.(eh1[2:end]./eh1[1:end-1])
    lgl2 =  [nothing; lgl2]
    lgh1 =  [nothing; lgh1]

    data = hcat(h, el2, eh1, lgl2, lgh1)
    header = ["h", L"$L_2$", L"$H_2$", L"$log_2(e_{L^2}) $", L"$log_2(e_{H_1}) $"]
    pretty_table(data, header=header, backend=Val(:latex ), formatters = ( ft_printf("%.3E"), ft_nonothing ))
end



