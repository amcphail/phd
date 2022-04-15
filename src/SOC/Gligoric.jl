using DifferentialEquations
using LinearAlgebra

using PyPlot

function unit_cell(u,u_m,u_p,B,λ,γ_m,γ_p,ζ)
    du = zeros(ComplexF64,length(u))

    du[1] = im*(B*u[1] + u[3] + u_p[3] + λ*(u_p[4] + im*u[4]) + (γ_p*abs2(u[1]) + ζ*abs2(u[2]))*u[1])
    du[2] = im*(-B*u[2] + u[4] + u_p[4] - λ*(u_p[3] - im*u[3]) + (γ_m*abs2(u[2]) + ζ*abs2(u[1]))*u[2])

    du[3] = im*(B*u[3] + u[1] + u_m[1] + u[5] + u_m[5] + λ*(u[6] - u_m[2] - im*(u[2] - u_m[6])) + (γ_p*abs2(u[3]) + ζ*abs2(u[4]))*u[3])
    du[4] = im*(-B*u[4] + u[2] + u_m[2] + u[6] + u_m[6] - λ*(u[5] - u_m[1] - im*(u[1] - u_m[5])) + (γ_m*abs2(u[4]) + ζ*abs2(u[3]))*u[4])

    du[5] = im*(B*u[5] + u[3] + u_p[3] - λ*(u[4] + im*u_p[4]) + (γ_p*abs2(u[5]) + ζ*abs2(u[6]))*u[5])
    du[6] = im*(-B*u[6] + u[4] + u_p[4] + λ*(u[3] - im*u_p[3]) + (γ_m*abs2(u[6]) + ζ*abs2(u[5]))*u[6])

    return du
end

function normalise(u,t,integrator)
    normalize!(u)

    return nothing
end

function normaliseCallback()
    return FunctionCallingCallback(normalise,func_everystep=true,func_start=true)
end

function diamond!(du,u,p,t)
    c = Int64(p[1])
    b = p[2]

    for i in 1:c
        v = u[((i-1)*6+1):i*6]
        if i > 1
            v_m = u[((i-2)*6+1):(i-1)*6]
        else
            if b == 0
                v_m = zeros(6)
            else
                v_m = u[((c-1)*6+1):c*6]
            end
        end
        if i < c
            v_p = u[(i*6+1):(i+1)*6]
        else
            if b == 0
                v_p = zeros(6)
            else
                v_p = u[1:6]
            end
        end

        d = unit_cell(v,v_m,v_p,p[3],p[4],p[5],p[6],p[7])

        for j in 1:6
            du[(i-1)*6+j] = d[j]
        end
    end

    return nothing
end

function latticeSolve(p)
    c = Int64(p[1])
    tspan = (p[8],p[9])

    u0 = zeros(ComplexF64,6*c)
    for i = 1:c
        u0[(i-1)*6+3] = 1
        u0[(i-1)*6+4] = 1
    end

    cb = normaliseCallback()

    #define the problem type
    prob = ODEProblem(diamond!,u0,tspan,p,callback=cb)

    sol = solve(prob,dt=0.000001)

    return sol
end

function runGligoric(b,l)
    p = [20, # unit cells
     b, # boundary
     0.0,  # B
     l,  # λ
     1.0,  # γ_p
     1.0,  # γ_m
     1.0,  # ζ
     0.0,  # t start
     100.0]  # t stop

    sol = latticeSolve(p)

    results = reduce(hcat,sol.u)'

    return results
end

function plotGligoric(fn,results)
    figure(figsize=(8,8))
    bar(1:120,3*abs2.(results))
#    xlim((0,0.4))
    title("Density")
    xlabel("Node")
    ylabel("Density")
    savefig(string("gligoric-",fn,".png"))
end

function plotGligorics(l,z,p)
    figure(figsize=(16,16))
    subplot(2,2,1)
    bar(1:120,3*abs2.(z))
#    xlim((0,0.4))
    title(string("Density, Lambda = ",l))
    xlabel("Node")
    ylabel("Density")
    subplot(2,2,2)
    bar(1:120,3*abs2.(p))
#    xlim((0,0.4))
    title(string("Density, Lambda = ",l))
    xlabel("Node")
    ylabel("Density")
    subplot(2,2,3)
    bar(54:66,3*abs2.(z[54:66]))
#    xlim((0,0.4))
    title("Density")
    xlabel("Node")
    ylabel("Density")
    subplot(2,2,4)
    bar(54:66,3*abs2.(p[54:66]))
#    xlim((0,0.4))
    title("Density")
    xlabel("Node")
    ylabel("Density")

    savefig(string("gligorics.png"))
end

function plotUnitCell(l,d)
    figure(figsize=(8,8))
    plot(3*abs2.(d[10:100,54]),label="a-")
    plot(3*abs2.(d[10:100,55]),label="a+")
    plot(3*abs2.(d[10:100,56]),label="b-")
    plot(3*abs2.(d[10:100,57]),label="b+")
    plot(3*abs2.(d[10:100,58]),label="c-")
    plot(3*abs2.(d[10:100,59]),label="c+")
#    xlim((0,0.4))
    title(string("Unit cell 9, Density Fluctuations, Lambda = ",l))
    xlabel("Time (a.u.)")
    ylabel("Density")
    legend()
    savefig(string("gligoric-","node-54-ts-",l,".png"))
end


λ = 10

z = runGligoric(0,λ)
p = runGligoric(1,λ)
plotGligorics(λ,z[end,:],p[end,:])
plotUnitCell(λ,z)
