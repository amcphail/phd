
module NaKPump

export pumpSolve

using DifferentialEquations

include("ForceFunctions.jl")

function hammer!(u,t,integrator)
    if u[4] != 0
        u[3] = 1
        u[4] = 0
    end
    return Nothing
end

f(out,in,t) = (out .= rand())

function hammerNoise()
    UniformNoiseFunction(t0) = NoiseFunction(t0,f,noise_prototype=rand(4))

    noise = UniformNoiseFunction(0)

    return noise
end

function hammerCallback()
    return FunctionCallingCallback(hammer!,func_everystep=true,func_start=true)
end

function damped_forced_periodic!(du,u,p,t,W)
    du[1] = u[2]
    pump = ForceFunctions.F_pumping(u[1],p[1],p[2],p[3])
    damp = ForceFunctions.F_damping(u[2],p[4])
    force = ForceFunctions.F_forcing(u[1],u[3],W[3],p[5],p[6])
    du[2] = pump + damp + force
    du[3] = - p[7]
    if force != 0
        du[4] = 1#eps(0.0)
    else
        du[4] = 0
    end

    return nothing
end

function pumpSolve(p)

    tspan = (p[8],p[9])

    u0 = zeros(4)
    u0[1] = (rand() - 0.5)

    noise = hammerNoise()

    cb = hammerCallback()

    #define the problem type
    prob = RODEProblem(damped_forced_periodic!,u0,tspan,p,noise=noise,callback=cb)

    #@time sol = solve(prob,RandomEM(),dt=1/100,saveat=t,force_dtmin=true);
    sol = solve(prob,RandomEM(),dt=1);

    return sol
end

end
