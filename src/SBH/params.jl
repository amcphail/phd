
using JLD2, FileIO

include("Connections.jl")

t0 = 0
tend = 1
dt = 0.01

T = [t0:dt:tend;]*10

nodes = 8
particles = 8
angles = 8

c = Connections.ring(nodes)

U = 1

ϵ = U # * particles

J_c = 1
J = J_c * U

φ = π/2

α = 0
β = 0

params = Dict([("T",T)
        ,("nodes",nodes)
        ,("particles",particles)
        ,("angles",angles)
        ,("edges",c)
        ,("ϵ",ϵ)
        ,("U",U)
        ,("J",J)
        ,("φ",φ)
        ,("α",α)
        ,("β",β)])

for i = 1:10
    params["J"] = i*particles
    params["α"] = 0

    jldsave(string("params-",i,".jld2");params)

    params["α"] = 1

    jldsave(string("params-stim-",i,".jld2");params)

    params["β"] = 1

    jldsave(string("params-phase-",i,".jld2");params)
end
