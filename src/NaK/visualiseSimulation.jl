using FileIO, JLD2

include("Visualisations.jl")

#use latex fonts for axis labelling
data_dir = "C:\\Users\\vivia\\Documents\\desktop\\NaKPump\\Julia\\data"
scratch_dir = data_dir

#simulation_id = ARGS[1]
#file_name = ARGS[2]

simulation_id = "003"
file_name = "C:\\Users\\vivia\\Documents\\desktop\\NaKPump\\Julia\\data\\003\\results-1-1-5.jld2"

pref = length(scratch_dir) + 1 + length(simulation_id) + 1

file_root = file_name[(begin+pref):(end-5)]

p_file = string(scratch_dir,"\\",simulation_id,"\\",file_root,".jld2")

jldopen(p_file,"r")
results = load(p_file,"results")
params = load(p_file,"params")
close(p_file)

vis_file = string(data_dir,"\\",simulation_id,"\\","visuals",file_root[(begin+7):end])

mi_fn = string(vis_file,".jld2")
png_fn = string(vis_file,".png")
gif_fn = string(vis_file,".gif")

jldopen(mi_fn,"r")
slice = load(mi_fn,"slice")
mi_f = load(mi_fn,"mi_f")
mi_b = load(mi_fn,"mi_b")
close(mi_fn)

Visualisations.plotNaKLatticeFour(png_fn,params[10],params[6],slice,mi_f,mi_b)
