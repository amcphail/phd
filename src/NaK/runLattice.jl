include("Connections.jl")
include("NaKLattice.jl")
include("Visualisations.jl")

n = 256

radius = 50
half_height = 15
pump_force = 0.0005
damping = 0.001
forcing_force = 0.01
forcing_prob = 0.1
relaxation = 0.3
lattice_spacing = 1 # 1
lattice_force = 0.0005

ti = 0.0
tf = 200000.0 #10000.0
tspan = (ti,tf)
step = 4
fraction = 20

connections = Connections.ring(n)

params = [n
        ,radius
        ,half_height
        ,pump_force
        ,damping
        ,forcing_force
        ,forcing_prob
        ,relaxation
        ,lattice_spacing
        ,lattice_force
        ,ti
		,tf
		,step
        ,fraction
		,connections]


@time sol = NaKLattice.latticeSolve(params)

results = reduce(hcat,sol.u)'

png_fn = "visuals.png"

Visualisations.prepareNaKLattice(png_fn,results,params)
#animateBars(gif_fn,slice)
