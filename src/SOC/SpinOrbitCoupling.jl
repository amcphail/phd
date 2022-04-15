
using QuantumOptics
using OrdinaryDiffEq, DiffEqCallbacks
using LinearAlgebra
using PyPlot

#use latex fonts for axis labelling
rc("text",usetex=true)
rc("font", family="serif")

ħ = 1.054571726e-34 # J ⋅ s
amu = 1.66054e-27 # kg
m = 86.9091835 * amu
k_b = 1.38064852e-23 # m^2 kg s^-2 K^-1

ω_z = 2*π*800.0 # rad s-1

a_0 = 5.29e-11 # m
a_s = 95*a_0
a_z = sqrt(ħ/(m*ω_z))

x_0 = 1e-6  # metres
#t_0 = 1e-6  # seconds
E_0 = (ħ^2/(2*m))/(x_0^2)
t_0 = ħ/E_0

E_k = E_0/k_b

#g = 4*π*ħ^2*a_s/m
#g_2D = (2*sqrt(2*π)*ħ^2*a_s)/(m*a_z)

ω_z_dimless = ω_z/t_0
a_s_dimless = a_s/x_0
a_z_dimless = sqrt(1/ω_z_dimless)
g_2D = (2*sqrt(2*π)*a_s_dimless)/(a_z_dimless)

k_dimless = 2*π*ω_z_dimless

m_dimless = 1

ΩR = 10. # pump strength
δ = 1e-8 # energy difference between hyperfine states
#k = 2*π # Raman beam wavevector
k = 0.5 # Raman beam wavevector

Npoints_x = 100
Npoints_y = 100

xmin = -10
xmax = 10

ymin = -10
ymax = 10

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

spinor_Tqp = Tqp ⊕ Tqp
spinor_Tpq = Tpq ⊕ Tpq

σ_x_12 = one(b_position)
σ_x_21 = one(b_position)

σ_x = σ_x_12 ⊕ σ_x_21
σ_x_11 = getblock(σ_x,1,2)
σ_x_22 = getblock(σ_x,2,1)
setblock!(σ_x,σ_x_12,1,2)
setblock!(σ_x,σ_x_21,2,1)
setblock!(σ_x,σ_x_11,1,1)
setblock!(σ_x,σ_x_22,2,2)

σ_x = σ_x/2

σ_z_11 = one(b_position)
σ_z_22 = -one(b_position)

σ_z = (σ_z_11 ⊕ σ_z_22)/2

H_kinetic = op_px^2/(2*m_dimless) + op_py^2/(2*m_dimless)

Ω_R_fun(x,y) = ΩR/2 * exp.(1im * 2 * k * y)
Ω_R = potentialoperator(b_position,Ω_R_fun)

#interaction strength
γ = 0.8
g = 10.0
g_uu = g
g_dd = g

g_du = g * (1 - γ)/(1 + γ)
g_ud = g_du

σ_pz = (op_py ⊕ (-op_py))/2

σ_pz_I = one(b_momentum) ⊕ (one(b_momentum))

#H1 = ΩR/2*σ_x + δ/2*σ_z

H_up = δ/2 * one(b_position)
H_down = - δ/2 * one(b_position)
H1 = H_up ⊕ H_down # diagonal blocks

setblock!(H1,Ω_R,1,2)
setblock!(H1,dagger(Ω_R),2,1)

Hψ_up = one(b_position)
Hψ_down = one(b_position)
Hψ = Hψ_up ⊕ Hψ_down

α = 0.5

trap(x,y) = α*(x^2 + y^2)

H_t = potentialoperator(b_position,trap)
H_trap = H_t ⊕ H_t

V_plot = Array(transpose(reshape(real.(diag(H_t.data)), (Npoints_x, Npoints_y))))

#spinor_H_kinetic = (H_kinetic ⊕ H_kinetic) -1im*(1/m_dimless)*k*σ_pz + (1/(2*m_dimless))*k^2*σ_pz_I
spinor_H_kinetic = (H_kinetic ⊕ H_kinetic) -1im*(1/m_dimless)*k*σ_pz + (1/(2*m_dimless))*k^2*σ_pz_I
spinor_H_kinetic_FFT = LazyProduct(spinor_Tqp,spinor_H_kinetic,spinor_Tpq)

H_total = LazySum(H1, H_trap,spinor_H_kinetic_FFT, Hψ)

# renormalization callback
norm_func(u, t, integrator) = normalize!(u)
ncb = FunctionCallingCallback(norm_func; func_everystep = true, func_start = true)

x0 = 0
y0 = 0
p0_x = 0
p0_y = 0
σ = 1

ψx0 = gaussianstate(b_position_x, x0, p0_x, σ)
ψy0 = gaussianstate(b_position_y, y0, p0_y, σ)

ψ0 = normalize!(ψx0 ⊗ ψy0)

spinor_ψ0 = normalize!(ψ0 ⊕ ψ0)

# integrator tolerances
abstol_int = 1e-8
reltol_int = 1e-8
maxiters_int = 1e8

# integration time settings
#t0 = 0
#tend = 20000
#dt = 1000;
t0 = 0
tend = 20
dt = 1;

dx = (xmax-xmin)/Npoints_x
dy = (ymax-ymin)/Npoints_y

steps = Npoints_x * Npoints_y

T = [t0:dt:tend;]*0.05

function H_im(t, ψ) # Update state-dependent term in H
    Hψ.data.nzval[1:steps] .= g_uu/(dx*dy) * abs2.(getblock(ψ,1).data) + g_ud/(dx*dy) * abs2.(getblock(ψ,2).data)
    Hψ.data.nzval[steps+1:end] .= g_dd/(dx*dy) * abs2.(getblock(ψ,2).data) + g_du/(dx*dy) * abs2.(getblock(ψ,1).data)
    return H_total
end

tout, ψt = timeevolution.schroedinger_dynamic(T, spinor_ψ0, H_im, callback=ncb)

density_up = [Array(transpose(reshape((abs2.(getblock(ψ,1).data)), (Npoints_x, Npoints_y)))) for ψ = ψt]
density_down = [Array(transpose(reshape((abs2.(getblock(ψ,2).data)), (Npoints_x, Npoints_y)))) for ψ = ψt]

xsample, ysample = samplepoints(b_position_x), samplepoints(b_position_y)

figure(figsize=(8, 8))
subplot(1,3,1)
title(L"\mathrm{spin}\uparrow, \, \, r-\mathrm{space}")
xlabel(L"x")
ylabel(L"y")
contourf(xsample, ysample, density_up[tend], cmap="Reds", alpha=1)
contourf(xsample, ysample, V_plot, alpha=0.4, cmap="Greys")

subplot(1,3,2)
title(L"\mathrm{spin}\downarrow, \, \, r-\mathrm{space}")
xlabel(L"x")
ylabel(L"y")
contourf(xsample, ysample, density_down[tend], cmap="Blues", alpha=1)
contourf(xsample, ysample, V_plot, alpha=0.4, cmap="Greys")

subplot(1,3,3)
title(L"\mathrm{spin}\uparrow,\downarrow, \, \, r-\mathrm{space}")
xlabel(L"x")
ylabel(L"y")
contourf(xsample, ysample, density_up[tend], cmap="Reds", alpha=0.5)
contourf(xsample, ysample, density_down[tend], cmap="Blues", alpha=0.5)
contourf(xsample, ysample, V_plot, alpha=0.4, cmap="Greys")

#tight_layout()


savefig("plot.pdf")

PyPlot.display_figs()
