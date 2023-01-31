module Results

    using Test
    using Plots
    import Dates
    import CairoMakie
    using LaTeXStrings
    using Latexify
    using PrettyTables

    function generate_plot(;ns, el2s, eh1s, ehs_energy, order=order, dirname=dirname)
        filename = dirname*"/conv_order_"*string(order)
        hs = 1 .// ns
        fig = CairoMakie.Figure()
        ax = CairoMakie.Axis(fig[1, 1], yscale = log10, xscale= log2,
                             yminorticksvisible = true, yminorgridvisible = true, yminorticks = CairoMakie.IntervalsBetween(8),
                             xlabel = L"h/{L}", ylabel = "error norms")

        CairoMakie.lines!(hs, el2s, label= L"$L^2$ norm", linewidth=2)
        CairoMakie.lines!(hs, eh1s, label= L"$H^1$ norm", linewidth=2)
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

    function generate_table(; ns, el2s, eh1s, ehs_energy, order=order, dirname=dirname)
        filename = dirname*"/conv_order_"*string(order)
        hs = 1 .// ns
        hs_str =  latexify.(hs)
        eoc_l2 = compute_eoc(hs, el2s)
        eoc_eh1 = compute_eoc(hs,eh1s)
        eoc_eh_energy = compute_eoc(hs,ehs_energy)

        println("==============")
        println("Order = $order")
        println("Mesh sizes = $hs")
        println("L2 errors  = $el2s")
        println("H1 errors  = $eh1s")
        println("Energy errors  = $eh1s")
        println("EOC L2 = $eoc_l2")
        println("EOC H1 = $eoc_eh1")
        println("EOC Energy = $eoc_eh_energy")

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

    function generate_figures(;ns, el2s, eh1s, ehs_energy, order::Integer, dirname::String)
        generate_plot(ns=ns, el2s=el2s, eh1s=eh1s, ehs_energy=ehs_energy, order=order, dirname=dirname)
        generate_table(ns=ns, el2s=el2s, eh1s=eh1s, ehs_energy=ehs_energy, order=order, dirname=dirname)
    end

end


