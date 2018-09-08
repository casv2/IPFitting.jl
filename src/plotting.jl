module Plotting

using Reexport
using StaticArrays, JuLIP, Plots, NBodyIPs

export slice_plot, force_plot, energy_plot

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

function slice_plot(IP)
    if Dict(IP)["components"][4]["D"]["__id__"] == "BondAngleDesc"
        E0 = IP.components[1].E0

        IP_comp = NBodyIPs.bodyorder.(IP.components)[2:end] #dropping E0 for 1body

        str = " body/cutoff: "
        for i in IP_comp
           str *= (string(i) * "/" * string(IP.components[i].D.cutoff.rcut)[1:4] * ", ")
        end

        #V2/V3/V4
        V2(r) = IP.components[2](r)
        V3(r1,r2,θ) = IP.components[3]( (SVector(r1, r2), SVector(cos(θ)*r1*r2) ) )
        #V4(r1, r2, θ) = IP.components[4]( (SVector(r1, r2), SVector(cos(θ)*r1*r2) ) )

        θ0 = acos(-1/3)
        rr2 = linspace(1, 8, 200)
        rr3 = linspace(1, 8, 200)

        P1 = plot( yaxis = ([-2, 0.7], ) )
        plot!(P1, rr2, V2.(rr2), label = "V2")
        plot!(P1, rr3, V3.(rr3, rr3, θ0), label = "V3(r,r,t0)")
        title!("e0=" * string(E0)[1:6] * str[1:end-2] * " BA")
        xlabel!("Interatomic distance (A)")
        ylabel!("Energy")
    end
end

function force_plot(test_data, IP; s = 50)
    data = sortcols(hcat([[split(test_data[i].configtype, ":")[1], test_data[i].D["F"], unfold(forces(IP, Atoms(test_data[i])))] for i in  1:s:length(test_data)]...))

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

    p1 = plot!(p1, [plot_min, plot_max], [plot_min, plot_max], color="black")

    plot(p1, p2, layout=(1,2), size=(1400,700))
end

function energy_plot(test_data, IP; s = 50)
    data = sortcols(hcat([[split(test_data[i].configtype, ":")[1], (test_data[i].D["E"][1]/length(test_data[i].at.Z)), (energy(IP, Atoms(test_data[i]))/length(test_data[i].at.Z))] for i in  1:s:length(test_data)]...))

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

    plot(p1, p2, layout=(1,2), size=(1400,700))
end

end
