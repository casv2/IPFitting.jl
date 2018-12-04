module Plotting
using Reexport

import JuLIP, NBodyIPFitting
using Plots, NBodyIPs

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
    IPs = NBodyIPs.bodyorder.(IP.components)[2:end]

    rr = linspace(xlim[1], xlim[2], 200)
    θ0 = acos(-1/3)

    p = plot( yaxis=([ylim[1],ylim[2]]) )

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

    if typeof(r0) == Float64
        vline!([r0], label="r0", color="black")
    end

    xlabel!("Interatomic distance (A)")
    ylabel!("Energy (eV)")

    if return_plot == false
        savefig("IP_plot.pdf")
    else
        display(p)
    end
end

end # module
