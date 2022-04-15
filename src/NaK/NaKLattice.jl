module NaKLattice

export latticeSolve

using DifferentialEquations

include("ForceFunctions.jl")

function hammer!(u,t,integrator)
    n = length(u)
    n = Int64(n/4)
    for i = 1:n
        if u[3*n+i] != 0
            u[2*n+i] = 1
            u[3*n+1] = 0
        end
    end
    return nothing
end

f(out,in,t) = (out .= rand())

function hammerNoise(n)
    UniformNoiseFunction(t0) = NoiseFunction(t0,f,noise_prototype=rand(4*n))

    noise = UniformNoiseFunction(0)

    return noise
end

function hammerCallback()
    return FunctionCallingCallback(hammer!,func_everystep=true,func_start=true)
end

function damped_forced_periodic!(du,u,p,t,W)
    n = Int64(p[1])
    for i = 1:n
        du[i] = u[n+i]
        pump = ForceFunctions.F_pumping(u[i],p[2],p[3],p[4])
        damp = ForceFunctions.F_damping(u[n+i],p[5])
        force = ForceFunctions.F_forcing(u[i],u[2*n+i],W[2*n+i],p[6],p[7])
        lattice = 0
        for j = 1:n
            c = p[15][i,j]
            if i != j  && c != 0
                lattice += c*ForceFunctions.F_lattice(u[i],u[j],p[9],p[10])
            end
        end
        du[n+i] = pump + damp + force + lattice
        du[2*n+i] = - p[8]
        if force != 0
            du[3*n+i] = 1
        else
            du[3*n+i] = 0
        end
    end
    return nothing
end

function latticeSolve(p)
    n = Int64(p[1])
    tspan = (p[11],p[12])

    u0 = zeros(4*n)
    for i = 1:n
        u0[i] = (rand() - 0.5)
    end

    noise = hammerNoise(n)

    cb = hammerCallback()

    #define the problem type
    prob = RODEProblem(damped_forced_periodic!,u0,tspan,p,noise=noise,callback=cb)

    len = (p[12] - p[11])/p[14]
    samples = Int64(len/p[13])

    front = LinRange(p[11],len,samples) |> collect
    back = LinRange(p[12]-len,p[12],samples) |> collect

    times = [front back]

    sol = solve(prob,RandomEM(),dt=1,saveat=times)

    return sol
end

end
