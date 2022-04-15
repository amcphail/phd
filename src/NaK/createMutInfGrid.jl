using Statistics

using FileIO, JLD2

include("settings.jl")

simulation_id = "002"

num_simulations = 20

file_root = string(data_dir,"/",simulation_id,"/")

cd(file_root)

mutinf = 0
sd = 0

mutinf_grid = zeros(num_simulations,num_simulations)
sd_grid = zeros(num_simulations,num_simulations)

x_b = 0
x_t = 0
y_b = 0
y_t = 0

for i in 1:num_simulations
    for j in 1:num_simulations
        global mutinf = 0
        global sd = 0
        for k in 1:5
            fn = string(file_root,"visuals-",i,"-",j,"-",k,".jld2")
            jldopen(fn)
            data = load(fn,"data")
            params = load(fn,"params")
            #close(fn)
            if i == 1
                global x_b = params[10]
            elseif i == num_simulations
                global x_t = params[10]
            end
            if j == 1
                global y_b = params[5]
            elseif j == num_simulations
                global y_t = params[5]
            end
            slice = data[end,1:100]
            s = std(slice)
            if (s > 0.025) && (maximum(slice) < 15) && (minimum(slice) > - 15)
                tot_b = load(fn,"tot_mi_b")
                global sd += s
            else
                tot_b = 0
            end
            global mutinf += tot_b
        end
        mutinf_grid[i,j] = mutinf / 5
        sd_grid[i,j] = sd / 5
    end
end

extent = [x_b,x_t,y_b,y_t]

jldsave(string(file_root,"mutinf.jld2"); mutinf_grid, sd_grid, extent)

