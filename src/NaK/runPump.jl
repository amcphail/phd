include("NaKPump.jl")
include("Visualisations.jl")

radius = 50
half_height = 20
pump_force = 0.0005
damping = 0.001
forcing_force = 0.01
forcing_prob = 0.1
relaxation = 0.001

ti = 0.0
tf = 50000.0 #10000.0
tspan = (ti,tf)

params = [radius
        ,half_height
        ,pump_force
        ,damping
        ,forcing_force
        ,forcing_prob
        ,relaxation
        ,ti
        ,tf]

@time sol = NaKPump.pumpSolve(params)

results = reduce(hcat,sol.u)'

png_fn = "pump.png"

Visualisations.plotNaKPump(png_fn,results,params)
