using DSP, Statistics

using PyPlot

include("MutualInformation.jl")

#use latex fonts for axis labelling
rc("text",usetex=true)
rc("font", family="serif")

linspace(a,b,n) = LinRange(a,b,n) |> collect

# data
#  psd
#  mi
#  animation
#  imshow
