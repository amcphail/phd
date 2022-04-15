using PyPlot
using Plots

using FileIO, JLD2

using LaTeXStrings, Revise

data_dir = "data"
simulation_id = "002"

file_root = string(data_dir,"/",simulation_id,"/")

cd(file_root)

fn = "mutinf.jld2"

jldopen(fn)
mutinf_grid = load(fn,"mutinf_grid")
extent = load(fn,"extent")
close(fn)

figure(figsize=(16,16))
imshow(mutinf_grid,extent=extent,vmin=0,vmax=7,aspect="auto")
title("Mutual Information")
xlabel("Lattice Force")
ylabel("Damping")
colorbar()


PyPlot.savefig(string(file_root,"mutinf.png"))
