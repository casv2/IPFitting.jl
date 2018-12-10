module Plotting
using Reexport

import JuLIP, NBodyIPFitting
using Plots, NBodyIPs, StaticArrays

gr(size=(800,500), html_output_format=:png)

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

function IP_plot(IP::NBodyIPs.NBodyIP; ylim = [-2,2], xlim = [0,8], r0 = 0, return_plot = true)
    # collect the IPs
    IPs = NBodyIPs.bodyorder.(IP.components)[2:end]

    rr = linspace(xlim[1], xlim[2], 200)
    θ0 = acos(-1/3)

    p = plot( yaxis=([ylim[1],ylim[2]]) )

    j = 2

    # plot V2a + V2b or V2

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

    # plot 3/4 BA/BL

    for i in deleteat!(IPs, findin(IPs, [2]))
        if i == 3
            try
                V3(r1,r2,r3) = IP.components[j]( SVector(r1,r2,r3) )
                plot!(p, rr, V3.(rr,rr,rr), label="V3 (BL)")
            catch
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

    # add r0

    if typeof(r0) == Float64
        vline!([r0], label="r0", color="black")
    end

    xlabel!("Interatomic distance (Angstrom)")
    ylabel!("Energy (eV)")

    if return_plot == false
        savefig("IP_plot.png")
    else
        display(p)
    end
end

function force_plot(IP::NBodyIPs.NBodyIP, test_data::Array{NBodyIPFitting.Dat,1}; s = 10, return_plot = true)
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
        savefig("forceplot.png")
    else
        display(P)
    end
end

function energy_plot(IP::NBodyIPs.NBodyIP, test_data::Array{NBodyIPFitting.Dat,1}; s = 10, return_plot = true)
    data = sortcols(hcat([[split(test_data[i].configtype, ":")[1], (test_data[i].D["E"][1]/length(test_data[i].at.Z)), (energy(IP, JuLIP.Atoms(test_data[i]))/length(test_data[i].at.Z))] for i in  1:s:length(test_data)]...))

    uniq_config_sel = unique(data[1,:])

    p1 = plot()
    p2 = plot()

    ylabel!(p1, "Predicted Energy per Atom (eV)")
    xlabel!(p1, "Target Energy per Atom(eV)")
    title!(p1, "Energy Scatter Plot")

    ylabel!(p2, "Predicted Energy per Atom (eV)")
    xlabel!(p2, "Target Energy per Atom (eV)")
    title!(p2, "Energy Error Plot")

    plot_max = findmax(unfold(data[3,:]))[1]
    plot_min = findmin(unfold(data[3,:]))[1]

    for uniq_config in uniq_config_sel
        indices = find(data[1,:] .== uniq_config)
        min, max = findmin(indices)[1], findmax(indices)[1]
        p1 = scatter!(p1, data[2, min:max], data[3, min:max], label=uniq_config, legend=:bottomright)
        p2 = scatter!(p2, data[2, min:max], abs.(data[3, min:max] - data[2, min:max]), label=uniq_config, legend=false, yaxis=(:log10, (0.00005,0.02)))
    end

    p1 = plot!(p1, [plot_min, plot_max], [plot_min, plot_max], color="black")

    P = plot(p1, p2, layout=(1,2), size=(1400,700))

    if return_plot == false
        savefig("energyplot.png")
    else
        display(P)
    end

end

function IP_pdf(IP::NBodyIPs.NBodyIP, info::Dict{String,Any})
    IP_plot(IP, return_plot = false)
    #error table
    error_table = "\\begin{supertabular}{ l c c c } \\toprule \n"
    error_table *= "Config type & E & F & V \\\\ \\midrule \n"

    for key in sort(collect(keys(info["errors"]["rmse"])))
        s = @sprintf "%s & %s & %s & %s \\\\ \n" replace(key, "_" => "\\_") string(info["errors"]["rmse"][key]["E"])[1:7] string(info["errors"]["rmse"][key]["F"])[1:7] string(info["errors"]["rmse"][key]["V"])[1:7]
        error_table *= s
    end
    error_table *= "\\end{supertabular}"

    #dataweights
    data_table = "\\begin{supertabular}{ l c c c } \\toprule \n"
    s = @sprintf "Data & %s & %s & %s \\\\ \\midrule \n" "E" "F" "V"
    data_table *= s
    s = @sprintf "Weight & %s & %s & %s \\\\ \\midrule \n" info["dataweights"]["E"] info["dataweights"]["F"] info["dataweights"]["V"]
    data_table *= s
    data_table *= "\\end{supertabular}"

    #weighttable
    weight_table = "\\begin{supertabular}{ l c } \\toprule \n"
    weight_table *= "Config type & Weight \\\\ \\midrule \n"

    for key in sort(collect(keys(info["configweights"])))
        s = @sprintf "%s & %s \\\\ \n" replace(key, "_" => "\\_") info["configweights"][key]
        weight_table *= s
    end
    weight_table *= "\\end{supertabular}"

    db = replace(info["dbpath"][3:end], "_" => "\\_")

    e0 = info["E0"]

    basis = string(length(info["Ibasis"]))

    template = "\\documentclass[a4paper,landscape]{article}
    \\usepackage{booktabs}
    \\usepackage[a4paper,margin=1in,landscape,twocolumn]{geometry}
    \\usepackage{amsmath}
    \\usepackage{graphicx}
    \\usepackage{subfig}
    \\usepackage{diagbox}
    \\usepackage{supertabular}
    \\begin{document}
    \\begin{center}
    \\textbf{LsqDB}: $db \\\\
    \\vspace{2mm}
    \\textbf{E0}: $e0 \\\\
    \\vspace{2mm}
    \\textbf{Basis functions}: $basis \\\\
    \\vspace{3mm}
    $data_table \\\\
    \\vspace{3mm}
    $weight_table \\\\
    \\vspace{3mm}
    $error_table \\\\
    \\begin{figure}[h]
        \\centering
        \\subfloat{{\\includegraphics[height=8cm]{IP_plot.png} }}%
        \\caption{Slices of \$V_{n}\$}%
    \\end{figure}
    \\end{center}
    \\end{document}"

    write("out.tex", template)

    run(`pdflatex out.tex`)
    run(`mv out.pdf PIP_analysis.pdf`)
    sleep(1)
    run(`rm out.tex out.log out.aux`)
end




end # module
