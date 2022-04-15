using Statistics
using StatsBase

using PyPlot

#use latex fonts for axis labelling
#rc("text",usetex=true)
#rc("font", family="serif")

function plotSpinorPhase(fn,V_plot,Npoints_x,Npoints_y,xsample,ysample,ψ)
    density_up = Array(transpose(reshape((abs2.(getblock(ψ,1).data)), (Npoints_x, Npoints_y))))
    density_down = Array(transpose(reshape((abs2.(getblock(ψ,2).data)), (Npoints_x, Npoints_y))))

    phase_up = Array(transpose(reshape((angle.(getblock(ψ,1).data)), (Npoints_x, Npoints_y))))
    phase_down = Array(transpose(reshape((angle.(getblock(ψ,2).data)), (Npoints_x, Npoints_y))))

    figure(figsize=(12, 12))
    subplot(2,2,1)
    #title("Spinor GPE with phase, spin up")
    title(L"\mathrm{Spinor\ GPE\ with\ Phase\ Shift} \, \, \mathrm{spin}\uparrow, \, \, |\Psi|^2")
    contourf(xsample, ysample, density_up, cmap="jet")
    contourf(xsample, ysample, V_plot, alpha=0.4, cmap="Greys")
    #annotate(xy=[25, 10], s="up", fontsize=20)

    subplot(2,2,2)
    #title("Spinor GPE with phase, spin down")
    title(L"\mathrm{Spinor\ GPE\ with\ Phase\ Shift} \, \, \mathrm{spin}\downarrow, \, \, |\Psi|^2")
    contourf(xsample, ysample, density_down, cmap="jet")
    contourf(xsample, ysample, V_plot, alpha=0.4, cmap="Greys")
    #annotate(xy=[25, 10], s="down", fontsize=20)

    subplot(2,2,3)
    #title(L"\mathrm{Spinor\ GPE} \, \, \mathrm{spin}\uparrow, \, \, \mathrm{phase}(\Psi)")
    imshow(phase_up,alpha=density_up ./ maximum(density_up), cmap="jet")
    #contourf(xsample, ysample, phase_up, cmap="jet")
    #contourf(xsample, ysample, V_plot, alpha=0.4, cmap="Greys")
    #annotate(xy=[25, 10], s="up", fontsize=20)

    subplot(2,2,4)
    #title(L"\mathrm{Spinor\ GPE} \, \, \mathrm{spin}\downarrow, \, \, \mathrm{phase}(\Psi)")
    imshow(phase_down,alpha=density_down ./ maximum(density_down), cmap="jet")
    #contourf(xsample, ysample, phase_down, cmap="jet")
    #contourf(xsample, ysample, V_plot, alpha=0.4, cmap="Greys")
    #annotate(xy=[25, 10], s="down", fontsize=20)


    savefig(fn)
end

function plotSpinor(fn,V_plot,Npoints_x,Npoints_y,xsample,ysample,ψ)
    density_up = Array(transpose(reshape((abs2.(getblock(ψ,1).data)), (Npoints_x, Npoints_y))))
    density_down = Array(transpose(reshape((abs2.(getblock(ψ,2).data)), (Npoints_x, Npoints_y))))

    phase_up = Array(transpose(reshape((angle.(getblock(ψ,1).data)), (Npoints_x, Npoints_y))))
    phase_down = Array(transpose(reshape((angle.(getblock(ψ,2).data)), (Npoints_x, Npoints_y))))

    figure(figsize=(12, 12))
    subplot(2,2,1)
    #title("Spinor GPE, spin up")
    title(L"\mathrm{Spinor\ GPE} \, \, \mathrm{spin}\uparrow, \, \, |\Psi|^2")
    contourf(xsample, ysample, density_up, cmap="jet")
    contourf(xsample, ysample, V_plot, alpha=0.4, cmap="Greys")
    #annotate(xy=[25, 10], s="up", fontsize=20)

    subplot(2,2,2)
    #title("Spinor GPE, spin down")
    title(L"\mathrm{Spinor\ GPE} \, \, \mathrm{spin}\downarrow, \, \, |\Psi|^2")
    contourf(xsample, ysample, density_down, cmap="jet")
    contourf(xsample, ysample, V_plot, alpha=0.4, cmap="Greys")
    #annotate(xy=[25, 10], s="down", fontsize=20)

    subplot(2,2,3)
    #title(L"\mathrm{Spinor\ GPE} \, \, \mathrm{spin}\uparrow, \, \, \mathrm{phase}(\Psi)")
    imshow(phase_up,alpha=density_up ./ maximum(density_up), cmap="jet")
    #contourf(xsample, ysample, phase_up, cmap="jet")
    #contourf(xsample, ysample, V_plot, alpha=0.4, cmap="Greys")
    #annotate(xy=[25, 10], s="up", fontsize=20)

    subplot(2,2,4)
    #title(L"\mathrm{Spinor\ GPE} \, \, \mathrm{spin}\downarrow, \, \, \mathrm{phase}(\Psi)")
    imshow(phase_down,alpha=density_down ./ maximum(density_down), cmap="jet")
    #contourf(xsample, ysample, phase_down, cmap="jet")
    #contourf(xsample, ysample, V_plot, alpha=0.4, cmap="Greys")
    #annotate(xy=[25, 10], s="down", fontsize=20)

    savefig(fn)
end

function plotPhases(fn,T,d)

    figure(figsize=(12,18))
    for j = 1:3

        N = length(d[j])

        μs_up = zeros(N)
        μs_down = zeros(N)

        σs_up = zeros(N)
        σs_down = zeros(N)

        ψt = d[j]

        for i = 1:N
            ψ = ψt[i]

            density_up = aweights(abs2.(getblock(ψ,1).data))
            density_down = aweights(abs2.(getblock(ψ,2).data))

            phase_up = angle.(getblock(ψ,1).data)
            phase_down = angle.(getblock(ψ,2).data)

            μs_up[i] = mean(phase_up,density_up)
            μs_down[i] = mean(phase_down,density_down)

            σs_up[i] = std(phase_up,density_up)
            σs_down[i] = std(phase_down,density_down)
        end

        subplot(3,2,j*2-1)
        title("Mean Phase")
        plot(T,μs_up) #,'r',label="up")
        #errorbars(1:N,σs_up)
        plot(T,μs_down) #,'b',label="down")
        #errorbars(1:N,σs_down)
        xlabel("Time")
        ylabel("Phase (radians)")
        legend()

        subplot(3,2,j*2)
        title("Std Dev. Phase")
        plot(T,σs_up) #,'r',label="up")
        plot(T,σs_down) #,'b',label="down")
        xlabel("Time")
        ylabel("Standard deviation")
        legend()

    end

    savefig(fn)
end

function plotDiamond(fn,V_plot,index)
    jldopen(string("diamond.jld2"),"r")

    xsample = load("diamond.jld2","xsample")
    ysample = load("diamond.jld2","ysample")
    Npoints_x = load("diamond.jld2","Npoints_x")
    Npoints_y = load("diamond.jld2","Npoints_y")
    ψt = load("diamond.jld2","ψt")

    density_up = [Array(transpose(reshape((abs2.(getblock(ψ,1).data)), (Npoints_x, Npoints_y)))) for ψ = ψt]
    density_down = [Array(transpose(reshape((abs2.(getblock(ψ,2).data)), (Npoints_x, Npoints_y)))) for ψ = ψt]

    figure(figsize=(8, 8))
    subplot(1,3,1)
    title("Spin Up")
    xlabel("x (um)")
    ylabel("y (um)")
#    title(L"\mathrm{spin}\uparrow, \, \, r-\mathrm{space}")
#    xlabel(L"x (\mu m)")
#    ylabel(L"y (\mu m)")
    contourf(xsample, ysample, density_up[index], cmap="Reds", alpha=1)
    contourf(xsample, ysample, V_plot, alpha=0.4, cmap="Greys")

    subplot(1,3,2)
    title("Spin Down")
    xlabel("x (um)")
#    title(L"\mathrm{spin}\downarrow, \, \, r-\mathrm{space}")
#    xlabel(L"x (\mu m)")
    #ylabel(L"y (\mu m)")
    contourf(xsample, ysample, density_down[index], cmap="Blues", alpha=1)
    contourf(xsample, ysample, V_plot, alpha=0.4, cmap="Greys")

    subplot(1,3,3)
    title(L"\mathrm{spin}\uparrow,\downarrow, \, \, r-\mathrm{space}")
    xlabel(L"x (\mu m)")
    #ylabel(L"y (\mu m)")
    contourf(xsample, ysample, density_up[index], cmap="Reds", alpha=0.5)
    contourf(xsample, ysample, density_down[index], cmap="Blues", alpha=0.5)
    contourf(xsample, ysample, V_plot, alpha=0.4, cmap="Greys")

    savefig(fn)

end
