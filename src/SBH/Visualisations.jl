using FileIO, JLD2

using PyPlot

jldopen("variances.jld2","r")
bh = load("variances.jld2","number")
close("variances.jld2")

jldopen("variances-stim.jld2","r")
sbh = load("variances-stim.jld2","number")
close("variances-stim.jld2")

jldopen("variances-phase.jld2","r")
pbh = load("variances-phase.jld2","number")
close("variances-phase.jld2")

ts = (0:100)./100

bhm = zeros(101,8,10)
sbhm = zeros(101,8,10)
pbhm = zeros(101,8,10)

for i in 1:10
    for j in 1:101
        run_bh = bh[i]
        run_sbh = sbh[i]
        run_pbh = pbh[i]
        bhm[j,:,i] = run_bh[j]
        sbhm[j,:,i] = run_sbh[j]
        pbhm[j,:,i] = run_pbh[j]
    end
end

index = 8
times = 1:5:101

figure(figsize=(24,8))
subplot(1,3,1)
imshow(bhm[times,:,index])
title("Bose Hubbard")
xlabel("Time (a.u.)")
ylabel("Variance of number over lattice sites")
subplot(1,3,2)
imshow(sbhm[times,:,index])
title("Stimulated Bose Hubbard")
xlabel("Time (a.u.)")
ylabel("Variance of number over lattice sites")
subplot(1,3,3)
imshow(pbhm[times,:,index])
title("Phase Bose Hubbard")
xlabel("Time (a.u.)")
ylabel("Variance of number over lattice sites")
#imshow(ss,vmin=0,vmax=6,extent=[0,4,0,1],origin='lower')
#ylabel("J")
#xlabel("ω")
#subplot(2,1,1)
#PyPlot.bar(1:b.N,ns)
#title(string("SBH: α: ",round(α,digits=3),", β: ",round(β,digits=3),", ω: ",round(ω,digits=3),", u: ",round(u,digits=3),", t: ",round(t,digits=3)))
#xlabel("Node")
#ylabel("Occupancy")
#subplot(2,1,2)
#PyPlot.bar(1:b.N,ts)
#title(string("SBH: α: ",round(α,digits=3),", β: ",round(β,digits=3),", ω: ",round(ω,digits=3),", U: ",round(U,digits=3),", J: ",round(J,digits=3)),", φ: ",round(φ,digits=3)))
#xlabel("Node")
#ylabel("Phase")

savefig("SBH3.png")

figure(figsize=(8,8))
imshow(sbhm[times,:,index],aspect="auto")
title("Stimulated Bose Hubbard, Number occupancy")
ylabel("Time (a.u.)")
xlabel("Lattice site")

savefig("SBH.png")

PyPlot.display_figs()
