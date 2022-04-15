
using QuantumOptics
using SparseArrays
using OrdinaryDiffEq, DiffEqCallbacks
using LinearAlgebra

using PyPlot

import LinearAlgebra: normalize!

include("Visualisations.jl")

ΩR = 0# 5. # pump strength
δ = 0#1e-10 # energy difference between hyperfine states
k = 2*π # Raman beam wavevector

Npoints_x = 100
Npoints_y = 100

xmin = -4
xmax = 4

ymin = -4
ymax = 4

b_position_x = PositionBasis(xmin,xmax,Npoints_x)
b_momentum_x = MomentumBasis(b_position_x)

b_position_y = PositionBasis(ymin,ymax,Npoints_y)
b_momentum_y = MomentumBasis(b_position_y)

b_position = b_position_x ⊗ b_position_y
b_momentum = b_momentum_x ⊗ b_momentum_y

op_x = position(b_position_x) ⊗ one(b_position_y)
op_y = one(b_position_x) ⊗ position(b_position_y)

op_px = momentum(b_momentum_x) ⊗ one(b_momentum_y)
op_py = one(b_momentum_x) ⊗ momentum(b_momentum_y)

Tqp = transform(b_position,b_momentum)
Tpq = transform(b_momentum,b_position)

H_kinetic = op_px^2/k^2 + op_py^2/k^2

Ω_R_fun(x,y) = ΩR/2 * exp.(1im * 2 * k * x)
Ω_R = potentialoperator(b_position,Ω_R_fun)

α = 1

trap_fun(x,y) = α*(2*x^2+y^2)
trap = potentialoperator(b_position,trap_fun)

spinor_trap = trap ⊕ trap

V_plot = Array(transpose(reshape(real.(diag(trap.data)), (Npoints_x, Npoints_y))))

xsample, ysample = samplepoints(b_position_x), samplepoints(b_position_y)

#interaction strength
γ = 0.8
g = 10.0
g_uu = g
g_dd = g

g_du = g * (1 - γ)/(1 + γ)
g_ud = g_du

spinor_Tqp = Tqp ⊕ Tqp
spinor_Tpq = Tpq ⊕ Tpq

spinor_H_kinetic = H_kinetic ⊕ H_kinetic
spinor_H_kinetic_FFT = LazyProduct(spinor_Tqp,spinor_H_kinetic,spinor_Tpq)

H_up = δ/2 * one(b_position)
H_down = - δ/2 * one(b_position)
H1 = H_up ⊕ H_down # diagonal blocks

setblock!(H1,Ω_R,1,2)
setblock!(H1,dagger(Ω_R),2,1)

Hψ_up = one(b_position)
Hψ_down = one(b_position)
Hψ = Hψ_up ⊕ Hψ_down

Hψx_up = one(b_position)
Hψx_down = one(b_position)
Hψx = Hψx_up ⊕ Hψx_down

Hϕ_up = one(b_position)
Hϕ_down = one(b_position)
Hϕ = Hϕ_up ⊕ Hϕ_down

H_U = LazySum(Hψ,LazyProduct(Hψx,Hϕ))

H_total = LazySum(H1, spinor_trap, spinor_H_kinetic_FFT) #, H_U)
H_gs = -1im*H_total

H = spinor_Tqp*spinor_H_kinetic*spinor_Tpq + H1 + spinor_trap
H = (H + dagger(H))/2

e_gs, ψ_gs = eigenstates(H,1,maxiter=1000000)

ψ_gs_up = getblock(ψ_gs[1],1)
ψ_gs_up.data = ψ_gs_up.data .* exp(im*rand()*2*π)

ψ_gs_down = getblock(ψ_gs[1],2)

spinor_ψ0 = ψ_gs_up ⊕ ψ_gs_down

# integrator tolerances
abstol_int = 1e-8
reltol_int = 1e-8
maxiters_int = 1e8

# integration time settings
t0 = 0
tend = 9.9
dt = 0.1;

dx = (xmax-xmin)/Npoints_x
dy = (ymax-ymin)/Npoints_y

steps = Npoints_x * Npoints_y

T = [t0:dt:tend;]*5

S = 1
ϕ = 0

function phase_shift_self(S,Npoints_x,Npoints_y,Ψ)
    data = Array(transpose(reshape(Ψ.data,(Npoints_x,Npoints_y))))
    phase = angle.(data)

    shift = zeros((Npoints_x,Npoints_y))

    for i in 1:Npoints_x
        for j in 1:Npoints_y
            if i > 1
                shift[i,j] += phase[i,j] - phase[i-1,j]
            end
            if i < Npoints_x
                shift[i,j] += phase[i,j] - phase[i+1,j]
            end
            if j > 1
                shift[i,j] += phase[i,j] - phase[i,j-1]
            end
            if j < Npoints_y
                shift[i,j] += phase[i,j] - phase[i,j+1]
            end
        end
    end
    θ = exp.(im .* (S .* cos.(shift .+ π/2) .* π))

    s = Vector(reshape(transpose(θ),(Npoints_x*Npoints_y)))

    d = spdiagm(0 => s)

    return SparseOperator(Ψ.basis,d)
end

function phase_shift_other(S,Npoints_x,Npoints_y,Ψ,Φ,ϕ)
    Ψ_data = Array(transpose(reshape(Ψ.data,(Npoints_x,Npoints_y))))
    Ψ_phase = angle.(Ψ_data)

    Φ_data = Array(transpose(reshape(Φ.data,(Npoints_x,Npoints_y))))
    Φ_phase = angle.(Φ_data)

    shift = (Ψ_phase - Φ_phase) ./ 2

    θ = exp.(im .* (S .* cos.(shift .+ ϕ) .* π))

    s = Vector(reshape(transpose(θ),(Npoints_x*Npoints_y)))

    d = spdiagm(0 => s)

    return SparseOperator(Ψ.basis,d)
end

function update_H(t, ψ)
    Hψ.data.nzval[1:steps] .= g_uu/(dx*dy) * abs2.(getblock(ψ,1).data)
    Hψ.data.nzval[steps+1:end] .= g_dd/(dx*dy) * abs2.(getblock(ψ,2).data)

    Hψx.data.nzval[1:steps] .= g_ud/(dx*dy) * abs2.(getblock(ψ,2).data)
    Hψx.data.nzval[steps+1:end] .= g_du/(dx*dy) * abs2.(getblock(ψ,1).data)

    θ_up = phase_shift_other(S,Npoints_x,Npoints_y,getblock(ψ,1),getblock(ψ,2),ϕ)
    θ_down = phase_shift_other(S,Npoints_x,Npoints_y,getblock(ψ,2),getblock(ψ,1),ϕ)

    Hϕ.data.nzval[1:steps] .= θ_up.data.nzval[1:steps]
    Hϕ.data.nzval[steps+1:end] .= θ_down.data.nzval[1:steps]
end

function H_gs_im(t, ψ) # Update state-dependent term in H
    update_H(t,ψ)

    return H_gs
end

function H_im(t, ψ) # Update state-dependent term in H
    update_H(t,ψ)

    return H_total
end

function normalize!(u, basis)
    ψ = Ket(basis,u)

    ψ_up = getblock(ψ,1)
    ψ_down = getblock(ψ,2)

    normalize!(ψ_up)
    normalize!(ψ_down)

    ψ = ψ_up ⊕ ψ_down

    for i = 1:length(ψ.data)
        u[i] = ψ.data[i]
    end
end

# renormalization callback
norm_func(u, t, integrator) = normalize!(u, basis(Hψ))

ncb = FunctionCallingCallback(norm_func; func_everystep = true, func_start = false)

@time tout, ψ_gs_t = timeevolution.schroedinger_dynamic([0:8:8;], spinor_ψ0, H_gs_im, callback=ncb)

ψ_it_up = getblock(ψ_gs_t[end],1)
ψ_it_up.data = ψ_it_up.data .* exp(im*rand()*2*π)

ψ_it_down = getblock(ψ_gs_t[end],2)

ψ0_it = ψ_it_up ⊕ ψ_it_down

@time tout, ψt_nog = timeevolution.schroedinger_dynamic(T, ψ0_it, H_im, callback=ncb)

H_U = LazySum(Hψ,Hψx)

H_total = LazySum(H1, spinor_trap, spinor_H_kinetic_FFT, H_U)

@time tout, ψt_spinor = timeevolution.schroedinger_dynamic(T, ψ0_it, H_im, callback=ncb)

H_U = LazySum(Hψ,LazyProduct(Hψx,Hϕ))

H_total = LazySum(H1, spinor_trap, spinor_H_kinetic_FFT, H_U)

@time tout, ψt_phase = timeevolution.schroedinger_dynamic(T, ψ0_it, H_im, callback=ncb)

for i in 1:length(T)
       plotSpinor(string("images/spinor-nog-",i,".png"),V_plot,Npoints_x,Npoints_y,xsample,ysample,ψt_nog[i])
       plotSpinor(string("images/spinor-",i,".png"),V_plot,Npoints_x,Npoints_y,xsample,ysample,ψt_spinor[i])
       plotSpinorPhase(string("images/spinor-phase-",i,".png"),V_plot,Npoints_x,Npoints_y,xsample,ysample,ψt_phase[i])
end

plotPhases(string("images/phase.png"),T,[ψt_nog,ψt_spinor,ψt_phase])
