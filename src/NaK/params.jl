using JLD2

include("settings.jl")

simulation_id = ARGS[1]

file_root = string(data_dir,"/",simulation_id,"/")
if !isdir(file_root)
        mkdir(file_root)
end

#num_simulations = 20

n = 100

radius = 50
half_height = 15
pump_force = 0.0005
damping = 0.001
forcing_force = 0.004
forcing_prob = 0.1
relaxation = 0.2
lattice_spacing = 1 # 1
lattice_force = 0.0003

ti = 0.0
tf = 200000.0 #10000.0
tspan = (ti,tf)
step = 5
fraction = 20

connections = zeros(n,n)
for i = 1:n
	for j = 1:n
		if abs(i-j) <= 1
			connections[i,j] = 1
		end
	end
end

#min_lattice_force = 0.0
#max_lattice_force = 0.01

#min_forcing_force = 0.0
#max_forcing_force = 0.002

#min_damping = 0
#max_damping = 0.004

#linspace(a,b,n) = LinRange(a,b,n) |> collect

#ps = LinRange(min_lattice_force,max_lattice_force,num_simulations) |> collect
#qs = LinRange(min_forcing_force,max_forcing_force,num_simulations) |> collect
#rs = LinRange(min_damping,max_damping,num_simulations) |> collect

#for i = 1:length(ps)
#    for j = 1:length(qs)
        params = [n
                ,radius
                ,half_height
                ,pump_force
                ,damping#rs[j]#damping
                ,forcing_force #qs[j] # forcing_force
                ,forcing_prob
                ,relaxation
                ,lattice_spacing
                ,lattice_force#ps[i] # lattice_force
                ,ti
                ,tf
		,step
                ,fraction
		,connections]

        fn = string(file_root,"parameters-1-1.jld2")
        jldsave(fn;fn,params)
#    end
#end
