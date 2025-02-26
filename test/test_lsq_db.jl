
using Test
using IPFitting, NBodyIPs, ProgressMeter, JuLIP
using NBodyIPs: blpolys
using JuLIP: decode_dict
Fit = IPFitting
DB = IPFitting.DB
Data = Fit.Data

function rand_data(sym, N, configtype="rand")
   cubic = (N == 1)
   at = bulk(sym, cubic=cubic) * N
   rattle!(at, 0.1)
   F = (N == 1) ? nothing : rand(3, length(at))
   return Dat(at, configtype; E = rand(), F = F, V = rand(3,3))
end

##
println("Double-Check (de-)dictionisation of basis: ")
basis1 = blpolys(2, ExpTransform(2.0, 3.0), CosCut(5.0, 7.0), 10)
basis2 = blpolys(3, ExpTransform(2.5, 3.0), CosCut2s(2.0, 2.5, 4.0, 5.5), 6)
println(@test decode_dict.( Dict.( basis1 ) ) == basis1)
println(@test decode_dict.( Dict.( basis2 ) ) == basis2)

println("Double-Check (de-)dictionisation of Dat: ")
data1 = [ rand_data(:Ti, 3, "md") for n = 1:10 ]
data2 = [ rand_data(:Ti, 1, "cell") for n = 1:10 ]
println(@test Dat.(Dict.(data1)) == data1)
println(@test Dat.(Dict.(data2)) == data2)

##
println("Create a temporary database.")
tmpdir = mktempdir()
dbpath = joinpath(tmpdir, "temp")
basis = [basis1; basis2]
data = [data1; data2]
db = nothing

try
   global db = DB.LsqDB(dbpath, basis, data)
   println(@test true)
catch
   @info("...something went wrong...")
   println(@test false)
end

println("checking consistency of db")
println(@test DB.dbpath(db) == dbpath)
println(@test DB.kronfile(dbpath) == dbpath * "_kron.h5")
println(@test DB.infofile(dbpath) == dbpath * "_info.json")
println(@test db.basis == basis)
println(@test db.configs == data)

println("re-load the database")
db1 = DB.LsqDB(dbpath)
println(@test DB.dbpath(db1) == dbpath)
println(@test db1.basis == db.basis)
println(@test db1.configs == db.configs)
println(@test db1.Ψ == db.Ψ)

##
println("re-load the database with mmap=true")
db2 = DB.LsqDB(dbpath, mmap=true)
println(@test DB.dbpath(db2) == dbpath)
println(@test db2.basis == db.basis)
println(@test db2.configs == db.configs)
println(@test db2.Ψ == db.Ψ)

##
println("Delete the temporary database")
rm(tmpdir; force=true, recursive=true)
