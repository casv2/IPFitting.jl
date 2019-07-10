
using JuLIP, Test, IPFitting, DataFrames
using JuLIP.Potentials: evaluate_d
using JuLIP.MLIPs: IPSuperBasis
using SHIPs
using IPFitting: Dat, LsqDB
Lsq = IPFitting.Lsq
Err = IPFitting.Errors
using LinearAlgebra: norm

# generate random data
function generate_data(species, L, rmax, N, calc)
   data = Dat[]
   for n = 1:N
      at = bulk(species; cubic=true, pbc=true) * L
      rattle!(at, rand() * rmax)
      E = energy(calc, at)
      F = forces(calc, at)
      push!(data, Dat(at, "rand"; E = E, F = F))
   end
   return data
end

r0 = rnn(:Si)
calc = StillingerWeber()
data = generate_data(:Si, 2, 0.33*r0, 100, calc)

# 2 stands for 2 neighbours i.e. body-order 3
b3basis(deg) = IPSuperBasis(
      PairBasis(deg, PolyTransform(2, r0), 2, cutoff(calc)),
      SHIPBasis(TotalDegree(deg, 1.5), 2, PolyTransform(3, r0), 2, 0.5*r0, cutoff(calc))
   )

##
err_erms = Float64[]
err_frms = Float64[]
degrees = [4, 8, 12, 16, 20]
for deg in degrees
   B = b3basis(deg)
   @show length(B)
   db = LsqDB("", B, data)
   IP, fitinfo = Lsq.lsqfit( db,
                         E0 = 0.0,
                         configweights = Dict("rand" => 1.0),
                         obsweights   = Dict("E" => 100.0, "F" => 1.0) )
   push!(err_erms, fitinfo["errors"]["relrmse"]["set"]["E"])
   push!(err_frms, fitinfo["errors"]["relrmse"]["set"]["F"])
end


##
df = DataFrame( :degrees => degrees,
                :relrms_E => err_erms,
                :relrms_F => err_frms )
display(df)

(@test minimum(err_erms) < 1e-4) |> println
(@test minimum(err_frms) < 1e-2) |> println
