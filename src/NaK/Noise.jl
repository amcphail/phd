module UniformNoise

using Statistics, DifferentialEquations

export UniformNoiseProcess, UniformNoiseProcess!, UniformNoiseFunction

@inline function UNIFORM_NOISE_DIST(W,dt,rng)
  if typeof(W.dw) <: AbstractArray && !(typeof(W.dW) <: SArray)
    return @fastmath rand(rng,W.dW)
  else
    return @fastmath rand(rng,typeof(W.dW))
  end
end

function INPLACE_UNIFORM_NOISE_DIST(rand_vec,W,dt,rng)
  rand!(rng,rand_vec)
end

function UNIFORM_NOISE_BRIDGE(W,W0,Wh,q,h,rng)
  if typeof(W.dW) <: AbstractArray
    return @fastmath q*abs(h)*rand(rng,W.dw)
  else
    return @fastmath q*abs(h)*rand(rng,typeof(W.dW))
  end
end

function INPLACE_UNIFORM_NOISE_BRIDGE(rand_vec,W,W0,Wh,q,h,rng)
  wiener_randn!(rng,rand_vec)
  #rand_vec .= sqrt((1.-q).*q.*abs(h)).*rand_vec.+q.*Wh
  sqrtcoeff = @fastmath sqrt((1-q)*q*abs(h))
  @. rand_vec = sqrtcoeff*rand_vec+q*Wh
end

UniformNoiseProcess(t0,W0,Z0=nothing) = NoiseProcess(t0,W0,Z0,UNIFORM_NOISE_DIST,UNIFORM_NOISE_BRIDGE)
UniformNoiseProcess!(t0,W0,Z0=nothing) = NoiseProcess(t0,W0,Z0,INPLACE_UNIFORM_NOISE_DIST,INPLACE_UNIFORM_NOISE_BRIDGE)

UniformNoiseFunction(t0) = NoiseFunction(t0,x -> rand,Z=nothing;reset=true)


end
