using DifferentialEquations

using FileIO, JLD2

include("NaKLattice.jl")

include("settings.jl")

simulation_id = ARGS[1]
repeats = parse(Int64,ARGS[2])
file_name = ARGS[3]

pref = length(data_dir) + 1 + length(simulation_id) + 1

file_root = file_name[(begin+pref):(end-5)]

p_file = string(data_dir,"/",simulation_id,"/",file_root,".jld2")

jldopen(p_file,"r")
params = load(p_file,"params")
#close(p_file)

# set parameter set
# Gibbs sampling
#  run ODE
#  slice of results
#  mutual information
#  provide error

for i in 1:repeats
        data_file = string("results",file_root[(begin+10):end],"-",i,".jld2")

        @time sol = NaKLattice.latticeSolve(params)

        results = reduce(hcat,sol.u)'

        jldsave(string(scratch_dir,"/",simulation_id,"/",data_file); data_file, params, results)
end
