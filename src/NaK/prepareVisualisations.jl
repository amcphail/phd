using Statistics
using DSP, FFTW

#using LaTeXStrings, Revise

using FileIO, JLD2

include("MutualInformation.jl")

#use latex fonts for axis labelling
#rc("text",usetex=true)
#rc("font", family="serif")

include("settings.jl")

simulation_id = ARGS[1]
file_name = ARGS[2]

#data_dir = "C:\\Users\\vivia\\Documents\\desktop\\NaKPump\\Julia\\data"
#scratch_dir = data_dir
#simulation_id = "003"
#file_name = "C:\\Users\\vivia\\Documents\\desktop\\NaKPump\\Julia\\data\\003\\results-1-1-1.jld2"

pref = length(scratch_dir) + 1 + length(simulation_id) + 1

file_root = file_name[(begin+pref):(end-5)]

p_file = string(scratch_dir,"/",simulation_id,"/",file_root,".jld2")

jldopen(p_file,"r")
results = load(p_file,"results")
params = load(p_file,"params")
#close(p_file)

(Nt,nodes) = size(results)

tim_1 = results[:,1]
mem_1 = results[:,3]
kck_1 = results[:,4]

len = Int64(floor(Nt/2))
n = Int64(nodes/4)

slice = results[(end-120+1):end,1:n]

data = results[(end-100+1):end,:]

front = Matrix{Float64}(undef,len,n)
back = Matrix{Float64}(undef,len,n)

psd_f = []
psd_b = []

frqs = []

for i in 1:n
    channel = results[:,i]

    fr = channel[begin:len]
    bk = channel[(end-len+1):end]

    front[:,i] = fr
    back[:,i] = bk

    pdg = Periodograms.welch_pgram(front[:,i])

    pw = Periodograms.power(pdg)

    push!(psd_f,pw)

    pdg = Periodograms.welch_pgram(back[:,i])

    pw = Periodograms.power(pdg)

    push!(psd_b,pw)

    push!(frqs,Periodograms.freq(pdg))
end

mi_f = MutualInformation.mutual_information_dataset(front,-10,10,100)
mi_b = MutualInformation.mutual_information_dataset(back,-10,10,100)

tot_mi_f = sum(mi_f)/(n^2)
tot_mi_b = sum(mi_b)/(n^2)

vis_file =string(data_dir,"/",simulation_id,"/","visuals",file_root[(begin+7):end])

mi_fn = string(vis_file,".jld2")

jldsave(mi_fn; params, data, slice, tot_mi_f, tot_mi_b, mi_f, mi_b, frqs, psd_f, psd_b)
