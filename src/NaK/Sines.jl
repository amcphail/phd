using Plots
using PyPlot

include("MutualInformation.jl")

#use latex fonts for axis labelling
rc("text",usetex=true)
rc("font", family="serif")

function calculateMI(k1=1,ω1=1,k2=0,ω2=0,nodes=100,steps=8000)
    dx = 1/nodes
    dt = 1/steps

    lattice = zeros(2,steps,nodes)

    for i = 1:steps
        for j = 1:nodes
            lattice[1,i,j] = sin(i*ω1*2*π*dt+j*k1*2*π*dx)*cos(i*ω2*2*π*dt+j*k2*2*π*dx)
            lattice[2,i,j] = k1*cos(i*ω2*2*π*dt+j*k2*2*π*dx) + (1-k1)*(2*rand()-1)
        end
    end

    mi1 = MutualInformation.mutual_information_dataset(lattice[1,:,:],-1,1,100)
    mi2 = MutualInformation.mutual_information_dataset(lattice[2,:,:],-1,1,100)

    return (mi1,mi2,lattice)
end

#(mi,lattice) = calculateMI(0,1,0,2)

ks = 0:0.1:1

len = length(ks)

ωs = 3*(1:len)/len

# mutinf = zeros(2,len,len)
#
# for i in 1:len
#     for j in 1:len
#
#         (mi1,mi2,lattice) = calculateMI(ks[i],1,0,ωs[j])
#
#         mutinf[1,i,j] = sum(mi1) ./ (100*8000)
#         mutinf[2,i,j] = sum(mi2) ./ (100*8000)
#     end
# end

subplot(1,2,1)
imshow(mutinf[1,:,:],extent=[0,1,3,0])
colorbar()
title(L"\mbox{Mutual Information (sines)}")
xlabel(L"k")
ylabel(L"\omega")
subplot(1,2,2)
imshow(mutinf[2,:,:],extent=[0,1,3,0])
colorbar()
title(L"\mbox{Mutual Information (with noise)}")
xlabel(L"\mbox{noise}")
ylabel(L"\omega")

PyPlot.savefig("mi-sines.png")
# skip = 100
# indices = 1:skip:steps
#
# slice = lattice[indices,:]
#
# subplot(1,2,1)
# imshow(slice')
# subplot(1,2,2)
# imshow(mi)

PyPlot.display_figs()
