"""
`module DB`

## Overview

This module implements a "database" for storing a precomputed LSQ system. This
is useful for fitting e.g. NBodyIPs since it allows one to precompute the <basis,
data> inner products which are very expensive, and then quickly construct LSQ
systems from them using e.g., many variants of weighting and regularisation.

A db called (e.g.) 'lsqdb' consists of two files:
 * `lsqdb_info.json` : this stores a list of basis functions and a list
 of configurations with attached observations (normally DFT data)
 * `lsqdb_kron.h5` : this stores all the inner products <basis, data>
 as a single large matrix. The observations and basis functions in
 `lsqdb_info` store the corresponding column and row indices.

## More details and remarks

### INFO

* BASIS: The basis is (for the time being) simply represented as a list (Vector) of
JuLIP.AbstractCalculator. These are stored as `INFO["basis"]`.

* DATA: Each "datum" is an atomistic configuration, which must have a
"configtype", stored as a `Dat`.
Informally two pieces of data with the same configtype should represent
"similar" kinds of configurations. `data = INFO["data"]` is a Dictionary
where the keys are the configtypes.

### (De-)Vectorising Data

Forces in `JuLIP` are represented as `Vector{SVector{...}}`, which is the
same memory layout as a 3 x N matrix (`JuLIP.mat`), which can can then be
vectorised (`[:]`), and this vectoriation is readily undone again.
Analogously, any data must be stored in such a vectorised format. This is
achieved (e.g. for forces) via
* `vec_obs(::Val{:F}, F) -> Vector`
* `devec_obs(::Val{:F}, Fvec) -> Vector{SVector{...}}`
or equivalently
* `vec_obs("F", F) -> Vector`
* `devec_obs("F", Fvec) -> Vector{SVector{...}}`
(the `Val` versions are used for performance optimisation)
"""
module DB

using Base.Threads:          SpinLock, nthreads
using StaticArrays:          SVector
using JuLIP:                 AbstractCalculator, AbstractAtoms, Atoms,
                             save_json, load_json, decode_dict
using IPFitting:        Dat, LsqDB, basis, eval_obs, observations,
                             observation, vec_obs, devec_obs,
                             tfor_observations
using IPFitting.Data:   configtype
using HDF5:                  h5open, read

import Base: flush, append!, union

export LsqDB, info, configtypes

const KRONFILE = "_kron.h5"
const INFOFILE = "_info.json"

"""
`dbpath(db::LsqDB)` : return the absolute path to the database files, not
including the 'kron.h5' and 'info.json' endings.
"""
dbpath(db::LsqDB) = db.dbpath

infofile(dbpath::AbstractString) = dbpath * INFOFILE

kronfile(dbpath::AbstractString) = dbpath * KRONFILE

# ------------ Save and load the info file

"if the file exists, append a random string to avoid overwriting"
function _backupfile(fname)
   if isfile(fname)
      fnew = fname * "." * String(rand('a':'z', 5))
      @warn("The file $fname already exists. It will be renamed to $fnew to avoid overwriting.")
      run(`mv $fname $fnew`)
   end
   return nothing
end

function load_info(dbpath::String)
   dbinfo = load_json(infofile(dbpath))
   basis = decode_dict.(dbinfo["basis"])
   configs = Dat.(dbinfo["configs"])   # here we already know the type
   return basis, configs
end

function save_info(dbpath::String, db)
   _backupfile(infofile(dbpath))
   save_json(infofile(dbpath),
             Dict("basis" => Dict.(db.basis),
                  "configs" => Dict.(db.configs))
            )
   return nothing
end

save_info(db) = save_info(dbpath(db), db)

# ----------- Save and load the KRON file

"save a single matrix to HDF5"
_savemath5(A, fname) =
   h5open(fname, "w") do fid
      fid["A"] = A
      nothing
   end

"load a single matrix from HDF5"
_loadmath5(fname) =
   h5open(fname, "r") do fid
      read(fid["A"])
   end

function save_kron(dbpath, db)
   _backupfile(kronfile(dbpath))
   _savemath5(db.Ψ, kronfile(dbpath))
end

save_kron(db) = save_kron(dbpath(db), db)

load_kron(dbpath::String; mmap=false) = _loadmath5(kronfile(dbpath))

load_kron(db::LsqDB; mmap=false) = load_kron(dbpath(db); mmap=mmap)

function LsqDB(dbpath::AbstractString; mmap=false)
   basis, configs = load_info(dbpath)
   Ψ = load_kron(dbpath; mmap=mmap)
   return LsqDB(basis, configs, Ψ, dbpath)
end


"""
for the time being, this just checks that the db director and db name are
admissible. In the future, this should initialise the DB files.
"""
function initdb(basedir, dbname)
   @assert !('/' in dbname)
   @assert isdir(basedir)
   dbpath = joinpath(basedir, dbname)
   initdb(dbpath)
   return nothing
end

function initdb(dbpath)
   # TODO: seems this fails to detect an existing database?
   @assert !isfile(infofile(dbpath))
   @assert !isfile(kronfile(dbpath))
   # check that we can actually create and delete this file
   save_json(infofile(dbpath), Dict("basis" => [], "data" => []))
   rm(infofile(dbpath))
   _savemath5(rand(5,5), kronfile(dbpath))
   rm(kronfile(dbpath))
   return nothing
end

function flush(db::LsqDB)
   save_info(db)
   save_kron(db)
   return nothing
end

function LsqDB(dbpath::AbstractString,
               basis::AbstractVector{<: AbstractCalculator},
               configs::AbstractVector{Dat};
               verbose=true,
               maxnthreads=nthreads())
   # assign indices, count observations and allocate a matrix
   Ψ = _alloc_lsq_matrix(configs, basis)
   # create the struct where everything is stored
   db = LsqDB(basis, configs, Ψ, dbpath)
   # parallel assembly of the LSQ matrix
   tfor_observations( configs,
      (n, okey, cfg, lck) -> safe_append!(db, lck, cfg, okey),
      msg = "Assemble LSQ blocks",
      verbose=verbose,
      maxnthreads=maxnthreads )
   # save to file
   if dbpath != ""
      verbose && @info("Writing db to disk...")
      try
         flush(db)
      catch
         @warn("""something went wrong trying to save the db to disk, but the data
               should be ok; if it is crucial to keep it, try to save manually.""")
      end
      verbose && @info("... done")
   else
      verbose && @info("db is not written to disk since `dbpath` is empty.")
   end
   return db
end


function set_matrows!(d::Dat, okey::String, irows::Vector{Int})
   d.rows[okey] = irows
   return d
end

matrows(d::Dat, okey::String) = d.rows[okey]

function _alloc_lsq_matrix(configs, basis)
   # the indices associated with the basis are simply the indices within
   # the array - there is nothing else to do here
   #
   # loop through all observations and assign indices
   nrows = 0
   for (okey, d, _) in observations(configs)
      len = length(observation(d, okey))
      set_matrows!(d, okey, collect(nrows .+ (1:len)))
      nrows += len
   end
   # allocate and return the matrix
   return zeros(Float64, nrows, length(basis))
end


function safe_append!(db::LsqDB, db_lock, cfg, okey)
   # computing the lsq blocks ("rows") can be done in parallel,
   lsqrow = evallsq(cfg, basis(db), okey)
   irows = matrows(cfg, okey)
   # but writing them into the DB must be done in a threadsafe way
   lock(db_lock)
   db.Ψ[irows, :] .= lsqrow
   unlock(db_lock)
   return nothing
end



# ------------------- Evaluating LSQ Blocks -----------------


"""
Take a basis and split it into individual basis groups.
"""
function split_basis(basis; splitfun = b -> hash(Val{:basis}(), b))
   # get some basis hashs of the individual basis functions
   tps = splitfun.(basis)
   Iord = Vector{Int}[]
   Bord = Any[]
   for tp in unique(tps)
      # find which elements of basis have type `tp`
      I = findall( [tp == t  for t in tps] )
      push!(Iord, I)
      push!(Bord, [b for b in basis[I]])
   end
   return Bord, Iord
end



# fill the LSQ system, i.e. evaluate basis at data points
function evallsq(d::Dat, B::AbstractVector{TB}, okey) where {TB <: AbstractCalculator}
   B1 = [b for b in B]
   if !(isconcretetype(eltype(B1)))
      return evallsq_split(d, B1, okey)
   end
   # TB is a leaf-type so we can use "evaluate_many"
   return _evallsq(Val(Symbol(okey)), B1, Atoms(d))
end

"""
evaluate one specific kind of datum, such as energy, forces, etc
for one Dat across a large section of the basis; implicitly this will
normally be done via `evaluate_many` or similar. The line
```
vals = eval_obs(vDT, B, at)
```
should likely be the bottleneck of `evallsq`
"""
function _evallsq(vDT::Val,
                  B::AbstractVector{<:AbstractCalculator},
                  at::AbstractAtoms)
   # vals will be a vector containing multiple evaluations
   vals = eval_obs(vDT, B, at)
   # vectorise the first so we know the length of the observation
   vec1 = vec_obs(vDT, vals[1])
   # create a multi-D array to reshape these into
   A = Array{Float64}(undef, length(vec1), length(B))
   A[:, 1] = vec1
   for n = 2:length(B)
      A[:, n] = vec_obs(vDT, vals[n])
   end
   return A
end

"""
split the Basis `B` into subsets of identical types and evaluate
those independently (fast).
"""
function evallsq_split(d, basis, okey)
   # TB is not a leaf-type so we should split the basis to be able to
   # evaluate_many & co
   Bord, Iord = split_basis(basis)
   lenobs = length(d.D[okey])
   lsqrow = zeros(lenobs, length(basis))
   for (B, IB) in zip(Bord, Iord)
      lsqrow[:, IB] = evallsq(d, B, okey)
   end
   return lsqrow
end


# =========================== Convenience =================================


# TODO: THIS NEEDS TO BE REWRITTEN! => ALSO WRITE TESTS!
# function union(db1::LsqDB,db2::LsqDB; dbpath = (db1.dbpath * "_u"))
#    @info("Warning: the union implies that the data in the two databases are the same")
#    configtypes = collect(keys(db1.data_groups))
#    basis = cat(1, db1.basis, db2.basis)
#    data_groups = db1.data_groups
#
#    kron_groups1 = db1.kron_groups
#    kron_groups2 = db2.kron_groups
#
#    kron_groups = Dict{String, KronGroup}()
#    for k in configtypes
#       kron_groups[k] = KronGroup()
#       @assert haskey(db2.data_groups, k)
#       for ot in ["E", "F", "V"]
#          if !haskey(data_groups[k][1].D, ot)
#             continue
#          end
#          kron_groups[k][ot] = cat(3,kron_groups1[k][ot],kron_groups2[k][ot])
#       end
#    end
#    db = LsqDB(basis, data_groups, kron_groups, dbpath)
#    flush(db)
#    return db
# end

# function union(db1::LsqDB,db2::LsqDB, dbpath = (db1.dbpath * "_u"))
#    basis = [basis(db1), basis(db2)]
#    data_groups::Dict{String, DataGroup}
#       kron_groups::Dict{String, KronGroup}
#       dbpath::String
#    end
#
#    basis(db::LsqDB) = db.basis
#    basis(db::LsqDB, i::Integer) = db.basis[i]
#    data(db::LsqDB) = db.data
# end


configtypes(db::LsqDB) = unique(configtype.(db.configs))

_nconfigs(db::LsqDB, ct::AbstractString) =
   length(find(configtype.(db.configs) .== ct))


function info(db::LsqDB)
   # config names, how many
   all_cts = configtype.(db.configs)
   unq_cts = unique(all_cts)
   configs_info = Dict{String, Int}()
   for ct in unq_cts
      configs_info[ct] = sum(all_cts .== ct)
   end

   println("======================================================")
   println("       LsqDB Summary ")
   println("------------------------------------------------------")
   println(" Datagroup: configname  =>  number of configs ")
   for (i, (key, val)) in enumerate(configs_info)
      println("        $i : $key  =>  $val ")
   end
   # basis groups, how many, max degree
   B = db.basis
   Bord, Iord = split_basis(B)
   println("------------------------------------------------------")
   println(" Basis Group  | length | description ")
   for (i, B1) in enumerate(Bord)
      lenstr = replace(string(length(B1), pad=5), '0' => ' ')
      idxstr = replace(string(i, pad=2), '0' => ' ')
      desc = string(B1)
      desc = desc[1:min(50, length(desc))]
      println("          $idxstr  | $lenstr  | $desc" )
      #
   end
   println("======================================================")
end

end
