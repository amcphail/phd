module Visualisations

export plotNaKPump
export prepareNaKLattice, plotNaKLattice
export plotNaKLatticeFour
export plotPhasePortrait

using DSP, FFTW
using Statistics

using LaTeXStrings

using PyPlot
#using Plots

include("MutualInformation.jl")
include("ForceFunctions.jl")

#use latex fonts for axis labelling
rc("text",usetex=true)
rc("font", family="serif")
rc("font", size=20)

linspace(a,b,n) = LinRange(a,b,n) |> collect

#function animateBars(fn,slice)
#    anim = @animate for i=1:Nt
#            Plots.bar(results[i,1:n],ylim=(-10,10))
#        end every skip
#    gif(anim,fn,fps=30)
#    return nothing
#end

function prepareNaKLattice(fn,results,params)

    x = params[10]
    y = params[6]

    (Nt,nodes) = size(results)

    len = Int64(floor(Nt/2))
    n = Int64(nodes/4)

    pos_1 = results[:,1]
    vel_1 = results[:,n+1]

    middle = results[Int64(floor((end-begin+1)/2)),1:n]
    last = results[end,1:n]

    slice = results[(end-256+1):end,1:n]

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

    position = results[:,1:n]
    velocity = results[:,(n+1):(2*n)]

    potential = position .^ 2
    kinetic = velocity .^ 2
    energy = potential .+ kinetic
    energy = abs.(energy .- mean(energy,dims=2))

    mi_f = MutualInformation.mutual_information_dataset(front,-15,15,150)
    mi_b = MutualInformation.mutual_information_dataset(back,-15,15,150)

    plotNaKLattice(fn,x,y,Nt,middle,last,slice,energy,pos_1,vel_1,frqs[1],sum(psd_f),sum(psd_b),mi_f,mi_b)
    #animateBars(gif_fn,slice)
end

function plotNaKLattice(fn,x,y,t,middle,last,slice,energy,pos,vel,freqs,pwr1,pwr2,mi1,mi2)
    (Nt,n) = size(slice)
    fig = figure(figsize=(24,36))
    subplot(3,3,1)
    imshow(slice',vmin=-5,vmax=5)
    title(string("Position of Oscillators, ","Lattice: ",round(x;digits=6),", Kick: ",round(y;digits=6)))
    xlabel("Time")
    ylabel("Node")
    subplot(3,3,9)
    PyPlot.loglog(freqs,pwr1,label="Initial")
    title("PSD")
    xlabel("Frequency")
    ylabel("Power")
    legend()
    subplot(3,3,9)
    PyPlot.loglog(freqs,pwr2,label="Final")
    title("PSD")
    xlabel("Frequency")
    ylabel("Power")
    legend()
    subplot(3,3,8)
    imshow(log.(abs2.(fft(slice)[Int64((n/2)+1):n,Int64((n/2)+1):n])))
    title("2D FFT")
    xlabel("Position")
    ylabel("Time")
#    PyPlot.plot(1:t,energy)
#    title("Energy of each oscillator")
#    xlabel("Time")
#    ylabel("Energy (std. dev.)")
    subplot(3,3,2)
    PyPlot.plot(1:t,pos)
    title("Position, Node 1")
    xlabel("Time")
    ylabel("Position")
    ylim((-15,15))
    subplot(3,3,5)
    PyPlot.plot(1:t,vel)
    title("Velocity, Node 1")
    xlabel("Time")
    ylabel("Velocity")
    ylim((-1,1))
    subplot(3,3,3)
    imshow(mi1,vmin=0,vmax=7)
    title("Mutual Information, beginning")
    xlabel("Node")
    ylabel("Node")
    subplot(3,3,6)
    imshow(mi2,vmin=0,vmax=7)
    title("Mutual Information, end")
    xlabel("Node")
    ylabel("Node")
    subplot(3,3,4)
    PyPlot.bar(1:n,middle)
    title(string("Equilibrating position, sigma=",round(std(middle),digits=3)," MI=",round(sum(mi1)/(100*100),digits=3)))
    xlabel("Node")
    ylabel("Position")
    ylim((-15,15))
    xlim((1,n))
    subplot(3,3,7)
    PyPlot.bar(1:n,last)
    title(string("Final position, sigma=",round(std(last),digits=3)," MI=",round(sum(mi2)/(100*100),digits=3)))
    xlabel("Node")
    ylabel("Position")
    ylim((-15,15))
    xlim((1,n))

    PyPlot.savefig(fn)

    return fig
end

function plotNaKPump(fn,results,params)
    (Nt,nodes) = size(results)

    ts = 1:Nt

    position = results[:,1]
    velocity = results[:,2]
    kick = results[:,4]

    fig = figure(figsize=(15, 12))
    subplot(2,2,1)
    PyPlot.plot(ts,position)
    title(L"Position")
    xlabel(L"time (t)")
    ylabel(L"position (z)")

    subplot(2,2,2)
    PyPlot.plot(ts,velocity)
    title(L"Velocity")
    xlabel(L"time (t)")
    ylabel(L"velocity (v)")

    subplot(2,2,3)
    PyPlot.plot(ts,kick)
    title(L"Kick")
    xlabel(L"time (t)")
    ylabel(L"kick (arb.)")

    subplot(2,2,4)
    zs = linspace(-20,20,1000)
    PyPlot.plot(zs,ForceFunctions.F_pumping.(zs,params[1],params[2],params[3]))

    PyPlot.savefig(fn)

    return fig
end

function plotNaKLatticeFour(fn,x,y,slice,mi1,mi2)
    fig = figure(figsize=(24,24))
    subplot(2,2,3)
    imshow(slice',vmin=-5,vmax=5,cmap="rainbow")
    title(string(L"Position of Oscillators, ",L"Lattice: ",round(x;digits=6),L", Damping: ",round(y;digits=6)))
    xlabel(L"Time")
    ylabel(L"Node")
    PyPlot.text(5,5,"(c)",color="k",fontsize=30)
    subplot(2,2,1)
    imshow(mi1,vmin=0,vmax=7)
    title(L"Mutual Information, beginning")
    xlabel(L"Node")
    ylabel(L"Node")
    PyPlot.text(5,5,"(a)",color="w",fontsize=30)
    subplot(2,2,2)
    imshow(mi2,vmin=0,vmax=7)
    title(L"Mutual Information, end")
    xlabel(L"Node")
    ylabel(L"Node")
    PyPlot.text(5,5,"(b)",color="w",fontsize=30)
    subplot(2,2,4)
    PyPlot.bar(1:100,slice[end,:],color="k")
    title(string(L"Final position, sigma=",round(std(slice[end,:]),digits=3),L",  MI=",round(sum(mi2)/(100*100),digits=3)))
    xlabel(L"Node")
    ylabel(L"Position")
    ylim((-5,5))
    xlim((1,100))
    PyPlot.text(5,4.5,"(d)",fontsize=30)

    PyPlot.savefig(fn)

    return fig
end

function plotPhasePortrait(fn,x,y)
    fig = figure(figsize=(8,8))
    title(L"Phase\ Portrait, Node\ 1")
    xlabel(L"Position")
    ylabel(L"Momentum")
    plot(x,y)

    PyPlot.savefig(fn)
end

function plotNode(fn,x,y)
    fig = figure(figsize=(8,8))
    title("Time Series, Node 1")
    xlabel("Time")
    ylabel("Position")
    plot(x,y)

    PyPlot.savefig(fn)
end
end
