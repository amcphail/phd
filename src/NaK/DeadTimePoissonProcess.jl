
using DifferentialEquations, Statistics, LaTeXStrings, Revise
using Random

using PyPlot
using Plots

include("Intervals.jl")

#use latex fonts for axis labelling
rc("text",usetex=true)
rc("font", family="serif")

linspace(a,b,n) = LinRange(a,b,n) |> collect

function F_forcing(armed,ran,k=1.,p=0.01)
    if armed < eps(0.0)
        if ran < p
            return - k
        end
    end
    return 0
end

function hammer!(u,t,integrator)
    if u[2] != 0
        u[1] = 1
        u[2] = 0
    end
    return Nothing
end

cb = FunctionCallingCallback(hammer!,func_everystep=true,func_start=false)

function dead_time!(du,u,p,t,W)
    force = F_forcing(u[1],W[1])
    du[1] = - p[1]
    if force != 0
        du[2] = 1
    else
        du[2] = 0
    end

    return nothing
end

u0 = zeros(2)

ti = 0.0
tf = 10000.0 #10000.0
tspan = (ti,tf)

decay = 0.1

params = [decay]

f(out,in,t) = (out .= rand())

UniformNoiseFunction(t0) = NoiseFunction(t0,f,noise_prototype=rand(2))
noise = UniformNoiseFunction(0)

#define the problem type
prob = RODEProblem(dead_time!,u0,tspan,params,noise=noise,callback=cb)

#@time sol = solve(prob,RandomEM(),dt=1/100,saveat=t,force_dtmin=true);
@time sol = solve(prob,RandomEM(),dt=1);

results = reduce(hcat,sol.u)'

(Nt,nodes) = size(results)

ts = 1:Nt

bins = 40

arm = results[:,1]
kick = results[:,2]

(e,w) = Intervals.create_hist(kick,bins)

figure(figsize=(15, 12))

subplot(3,1,1)
PyPlot.plot(ts,arm)
title(L"Arm")
xlabel(L"time (t)")
ylabel(L"arm (arb.)")

subplot(3,1,2)
PyPlot.plot(ts,kick)
title(L"Kick")
xlabel(L"time (t)")
ylabel(L"kick (arb.)")

subplot(3,1,3)
PyPlot.bar(e[1][begin:(end-1)],w,width=5)
title(L"Wait Times")
xlabel(L"Interval (t)")
ylabel(L"Count")

PyPlot.display_figs()
