module Plotting

using Reexport
using NBodyIPs, DataFrames, NBodyIPFitting, StaticArrays, JuLIP, Plots

export force_plot, energy_plot, IP_table, IP_info, IP_plot, IP_pdf

function unfold(A)
    V = []
    for x in A
        if x === A
            push!(V, x)
        else
            append!(V, unfold(x))
        end
    end
    V
end

function IP_table(info)
    df = DataFrame(config = [], E = [], F = [], V = [])

    for key1 in keys(info["errors"]["rmse"])
        df2 = DataFrame(config = [key1], E = [info["errors"]["rmse"][key1]["E"]],
                       F = [info["errors"]["rmse"][key1]["F"]], V = [info["errors"]["rmse"][key1]["V"]])
        append!(df, df2)
    end

    names!(df, [Symbol("config type"), Symbol("E (eV)"), Symbol("F (eV/A)"), Symbol("V (eV/A)")])

    rmse_table = reprmime("text/latex", df)

    df = DataFrame(config = [], E = [], F = [], V = [])

    for key1 in keys(info["errors"]["mae"])
        df2 = DataFrame(config = [key1], E = [info["errors"]["mae"][key1]["E"]],
                       F = [info["errors"]["mae"][key1]["F"]], V = [info["errors"]["mae"][key1]["V"]])
        append!(df, df2)
    end

    names!(df, [Symbol("config type"), Symbol("E (eV)"), Symbol("F (eV/A)"), Symbol("V (eV/A)")])

    mae_table = reprmime("text/latex", df)

    df_1 = DataFrame(config = [], weight = [])

    for key in keys(info["configweights"])
        df2 = DataFrame(config = key, weight = info["configweights"][key])
        append!(df_1, df2)
    end

    df_2 = DataFrame(data = [], weight = [])

    for key in keys(info["dataweights"])
        df2 = DataFrame(data = key, weight = info["dataweights"][key])
        append!(df_2, df2)
    end

    y = reprmime("text/latex", df_1)
    x = reprmime("text/latex", df_2)

    rmse = "RMSE"
    mae = "MAE"

    dbpath = replace(split(info["dbpath"], "/")[end], "_" => " ")

    E0 = info["E0"]

    if length(info["regularisers"]) == 0
        regulariser = "None"
    end

    b = "\\documentclass[preview]{standalone}
    \\usepackage{graphicx}
    \\usepackage{caption}
    \\captionsetup[figure]{font=small}
    \\begin{document}
    \\begin{table}
    \\caption{$rmse}
    \\begin{center}
    \\scalebox{1}{
    $rmse_table}
    \\end{center}
    \\end{table}
    \\begin{table}
    \\caption{$mae}
    \\begin{center}
    \\scalebox{1}{
    $mae_table}
    \\end{center}
    \\end{table}
    \\begin{center}
    \\begin{tabular}{ll}
    \\begin{tabular}{ccc}
    $x
    \\end{tabular}
    &
    \\begin{tabular}{ccc}
    $y
    \\end{tabular}
    \\end{tabular}
    \\newline
    \\newline
    database: $dbpath
    \\newline
    E0 = $E0
    \\newline
    regulariser : $regulariser
    \\newline
    \\end{center}
    \\end{document}"

    write("IP_table.tex", b)

    run(`pdflatex IP_table.tex`)
end

function IP_plot(IP; return_plot = true)
    IPs = NBodyIPs.bodyorder.(IP.components)[2:end]

    rr = linspace(1, 8, 200)
    θ0 = acos(-1/3)

    p = plot( yaxis=([-1,1]) )

    j = 2

    if length(find(IPs .== 2)) > 1
        V2a(r) = IP.components[j](r)
        j += 1
        V2b(r) = IP.components[j](r)
        plot!(p, rr, V2b.(rr) + V2a.(rr), label="V2a + V2b")
    else
        V2(r) = IP.components[j](r)
        plot!(p, rr, V2.(rr), label="V2")
    end

    j += 1

    for i in deleteat!(IPs, findin(IPs, [2]))
        println(i)
        if i == 3
            try
                println(j)
                V3(r1,r2,r3) = IP.components[j]( SVector(r1,r2,r3) )
                plot!(p, rr, V3.(rr,rr,rr), label="V3 (BL)")
            catch
                println(j)
                V3(r1,r2,θ) = IP.components[j]( (SVector(r1, r2), SVector(cos(θ)*r1*r2) ) )
                plot!(p, rr, V3.(rr, rr, θ0), label = "V3 (BA)")
            end
        elseif i == 4
            try
                V4(r1,r2,r3,r4,r5,r6) = IP.components[j]( SVector(r1,r2,r3,r4,r5,r6) )
                plot!(p, rr, V4.(rr,rr,rr,rr,rr,rr), label="V4 (BL)")
            catch
                V4(r1,r2,r3, θ1, θ2, θ3) = IP.components[j]( (SVector(r1, r2, r3), SVector(θ1, θ2, θ3)) )
                plot!(p, rr, V4.(rr, rr, rr, θ0, θ0, θ0), label="V4 (BA)")
            end
        end
        j += 1
    end

    r0 = rnn(:Si)

    vline!([r0], label="r0", color="black")
    xlabel!("Interatomic distance (A)")
    ylabel!("Energy (eV)")

    if return_plot == false
        savefig("IP_plot.pdf")
    else
        display(p)
    end
end

function IP_info(info; return_plot = true)
    b = bar(legend=:false)

    for key in keys(info["numconfigs"])
        bar!(b, [key], [info["numconfigs"][key]])
    end

    xlabel!("Configs in fit")
    ylabel!("Number of configs")

    if return_plot == false
        savefig("IP_info.pdf")
    else
        display(b)
    end
end

function force_plot(test_data, IP; s = 50, return_plot = true)
    data = sortcols(hcat([[split(test_data[i].configtype, ":")[1], test_data[i].D["F"], unfold(forces(IP, JuLIP.Atoms(test_data[i])))] for i in  1:s:length(test_data)]...))
    uniq_config_sel = unique(data[1,:])

    p1 = plot()
    p2 = plot()

    ylabel!(p1, "Predicted Force (eV/A)")
    xlabel!(p1, "Target Force (eV/A)")
    title!(p1, "Force Scatter Plot")

    ylabel!(p2, "Predicted Force Error (eV/A)")
    xlabel!(p2, "Target Force (eV/A)")
    title!(p2, "Force Error Plot")

    plot_max = findmax(unfold(data[3,:]))[1]
    plot_min = findmin(unfold(data[3,:]))[1]
    for uniq_config in uniq_config_sel
        indices = find(data[1,:] .== uniq_config)
        min, max = findmin(indices)[1], findmax(indices)[1]
        p1 = scatter!(p1, unfold(data[2, min:max]), unfold(data[3, min:max]), label=uniq_config, legend=:bottomright)
        p2 = scatter!(p2, unfold(data[2, min:max]), abs.(unfold(data[3, min:max]) - unfold(data[2, min:max])), label=uniq_config, legend=false, yaxis=(:log10, (0.00005,2)))
    end

    P = plot!(p1, [plot_min, plot_max], [plot_min, plot_max], color="black")

    P = plot(p1, p2, layout=(1,2), size=(1400,700))

    if return_plot == false
        savefig("IP_force.pdf")
    else
        display(P)
    end
end

function energy_plot(test_data, IP; s = 50, return_plot = true)
    data = sortcols(hcat([[split(test_data[i].configtype, ":")[1], (test_data[i].D["E"][1]/length(test_data[i].at.Z)), (energy(IP, JuLIP.Atoms(test_data[i]))/length(test_data[i].at.Z))] for i in  1:s:length(test_data)]...))

    uniq_config_sel = unique(data[1,:])

    p1 = plot()
    p2 = plot()

    ylabel!(p1, "Predicted Energy per Atom (eV)")
    xlabel!(p1, "Target Energy per Atom(eV)")
    title!(p1, "Energy Scatter Plot")

    ylabel!(p2, "Predicted Force Error (ev/A)")
    xlabel!(p2, "Target Energy per Atom (ev/A)")
    title!(p2, "Energy Error Plot")

    plot_max = findmax(unfold(data[3,:]))[1]
    plot_min = findmin(unfold(data[3,:]))[1]

    for uniq_config in uniq_config_sel
        indices = find(data[1,:] .== uniq_config)
        min, max = findmin(indices)[1], findmax(indices)[1]
        p1 = scatter!(p1, data[2, min:max], data[3, min:max], label=uniq_config, legend=:bottomright)
        p2 = scatter!(p2, data[2, min:max], abs.(data[3, min:max] - data[2, min:max]), label=uniq_config, legend=false, yaxis=(:log10, (0.00005,0.5)))
    end

    p1 = plot!(p1, [plot_min, plot_max], [plot_min, plot_max], color="black")

    P = plot(p1, p2, layout=(1,2), size=(1400,700))

    if return_plot == false
        savefig("IP_energy.pdf")
    else
        display(P)
    end

end

function IP_pdf(IP, info, test_data)
    IP_table(info)
    IP_plot(IP, return_plot = false)
    IP_info(info, return_plot = false)
    force_plot(test_data, IP, s=50, return_plot = false)
    energy_plot(test_data, IP, s=50, return_plot = false)

    #run(`convert IP_energy.png IP_energy.pdf`)
    #run(`convert IP_force.png IP_force.pdf`)
    run(`pdfjoin IP_plot.pdf IP_info.pdf IP_table.pdf IP_energy.pdf IP_force.pdf --outfile IP_pdf.pdf --rotateoversize false`)
end

end
