
using Random

using PyPlot

using JLD2, FileIO

include("StimulatedBoseHubbard.jl")

fn = ARGS[1]

param_dir = "src/SBH"

data_dir = "data/SBH"

ending = fn[(length(param_dir)+8):end]

Random.seed!(0)

jldopen(fn,"r")
params = load(fn,"params")

b = StimulatedBoseHubbard.LatticeBasis(params["nodes"],params["particles"],params["angles"])

Ψ0 = StimulatedBoseHubbard.rand(b,params["particles"])

@time tout, Ψt = StimulatedBoseHubbard.latticeSolve(params,Ψ0)

o_fn = string(data_dir,"/","results",ending)
jldsave(o_fn;params,tout,Ψt)

