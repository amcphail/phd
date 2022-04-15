
using QuantumOptics
using OrdinaryDiffEq, DiffEqCallbacks
using LinearAlgebra
using PyPlot

using JLD2, FileIO

include("Visualisations.jl")

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

ω_z_dimless = ω_z*t_0
a_s_dimless = a_s/x_0
a_z_dimless = sqrt(1/ω_z_dimless)
g_2D = (2*sqrt(2*π)*a_s_dimless)/(a_z_dimless)

#k_dimless = 2*π/ω_z_dimless
k_dimless = (2*π/(5.845e6*t_0))/x_0
E_z = 5.845e6*t_0

m_dimless = 1

E_Ω = E_k*k_dimless

Ω = 4*E_Ω#10. # pump strength
δ = E_z#1e-8 # energy difference between hyperfine states
#k = 2*π # Raman beam wavevector
#k = k_dimless # Raman beam wavevector

#interaction strength
γ = 0.8
g = 10.0
g = g_2D

g_uu = g
g_dd = g

g_du = g * (1 - γ)/(1 + γ)
g_du = g
g_ud = g_du

distance = 0.8
separation = 1.5
height = round(30e-9/E_k) # nK

Npoints_x = 100
Npoints_y = 100

xmin = -4
xmax = 4

ymin = -7
ymax = 7

# integrator tolerances
abstol_int = 1e-8
reltol_int = 1e-8
maxiters_int = 1e10

# integration time settings
#t0 = 0
#tend = 20000
#dt = 1000;
t0 = 0
tend = 0.1
dt = 0.01

dx = ((xmax-xmin)/Npoints_x)
dy = ((ymax-ymin)/Npoints_y)

steps = Npoints_x * Npoints_y

T = [t0:dt:tend;]

x0 = 0
y0 = 0
p0_x = 0
p0_y = 0
σ = 1

b_position_x = PositionBasis(xmin,xmax,Npoints_x)
b_momentum_x = MomentumBasis(b_position_x)

b_position_y = PositionBasis(ymin,ymax,Npoints_y)
b_momentum_y = MomentumBasis(b_position_y)

b_position = b_position_x ⊗ b_position_y
b_momentum = b_momentum_x ⊗ b_momentum_y

Hψ_up = one(b_position)
Hψ_down = one(b_position)
Hψ = Hψ_up ⊕ Hψ_down

function insideUnitCell(x,y,x0,d,s)
    if ((sqrt((x+x0+1+s)^2 + (y+s)^2) <= d)
        || (sqrt((x+x0+1+s)^2 + (y-s)^2) <= d)
        || (sqrt((x+x0+s)^2 + y^2) <= d))
       return true
   else
       return false
   end
end

function insideNE(x,y,x0,y0,s,d)
    if ((y - y0 < (x - x0) + d
        && (y - y0 > (x - x0) - d))
        && (y - y0 > - (x - x0))
        && (y - y0 < - (x - x0) + sqrt(2)*s))
        return true
    else
        return false
    end
end

function insideNW(x,y,x0,y0,s,d)
    y0 = -y0
    if ((y - y0 < - (x - x0) + 2*s + d
        && (y - y0 > - (x - x0) + 2*s - d))
        && (y - y0 > (x - x0))
        && (y - y0 < (x - x0) + sqrt(2)*s))
        return true
    else
        return false
    end
end

function insideUnitCellEdges(x,y,x0,d,s)
    if ((sqrt((x-x0+s+s)^2 + (y+s)^2) <= d)
        || insideNE(x,y,x0-2*s,-s,s,0.25)
        || (sqrt((x-x0+s+s)^2 + (y-s)^2) <= d)
        || insideNW(x,y,x0-2*s,s,s,0.25)
        || (sqrt((x-x0+s)^2 + y^2) <= d)
        || insideNE(x,y,x0-s,0,s,0.25)
        || insideNW(x,y,x0-s,2*s,s,0.25))
       return true
   else
       return false
   end
end

function wells(y, x, d, s, h)

#    if ((x > -4 && x < 0) || ((x > d) && (x < (4 + d)))) && (y > -2 && y < 2)
#        return 0
    if (insideUnitCellEdges(x,y,4.5,d,s)
        || insideUnitCellEdges(x,y,1.5,d,s)
        || insideUnitCellEdges(x,y,-1.5,d,s)
        || (sqrt((x-6+s)^2 + (y+s)^2) <= d)
        || (sqrt((x-6+s)^2 + (y-s)^2) <= d))
        return 0
    else
        return h
    end
end

potential(x,y) = wells(x,y,distance,separation,height)

V = potentialoperator(b_position, potential)

V_trap = V ⊕ V

V_plot = Array(transpose(reshape(real.(diag(V.data)), (Npoints_x, Npoints_y))))

function initialState(x_basis,y_basis, x0, y0, p0_x, p0_y, σ, s)

    ψx0 = gaussianstate(x_basis, x0, p0_x, σ)
    ψy0 = gaussianstate(y_basis, y0, p0_y, σ)

    ψx1 = gaussianstate(x_basis, x0, p0_x, σ)
    ψy1 = gaussianstate(y_basis, y0 + 2*s, p0_y, σ)

    ψx2 = gaussianstate(x_basis, x0, p0_x, σ)
    ψy2 = gaussianstate(y_basis, y0 - 2*s, p0_y, σ)

    ψ0 = normalize!(ψx0 ⊗ ψy0 + ψx1 ⊗ ψy1 + ψx2 ⊗ ψy2)

    spinor_ψ0 = normalize!(ψ0 ⊕ ψ0)

    return spinor_ψ0
end

function makeHamiltonian(qx,qy,q,px,py,p,m,k,Ω,δ,V,Hψ)
    op_x = position(qx) ⊗ one(qy)
    op_y = one(qx) ⊗ position(qy)

    op_px = momentum(px) ⊗ one(py)
    op_py = one(px) ⊗ momentum(py)

    Tqp = transform(q,p)
    Tpq = transform(p,q)

    spinor_Tqp = Tqp ⊕ Tqp
    spinor_Tpq = Tpq ⊕ Tpq

    σ_x_12 = one(q)
    σ_x_21 = one(q)

    σ_x = σ_x_12 ⊕ σ_x_21
    σ_x_11 = getblock(σ_x,1,2)
    σ_x_22 = getblock(σ_x,2,1)
    setblock!(σ_x,σ_x_12,1,2)
    setblock!(σ_x,σ_x_21,2,1)
    setblock!(σ_x,σ_x_11,1,1)
    setblock!(σ_x,σ_x_22,2,2)

    σ_x = σ_x/2

    σ_z_11 = one(q)
    σ_z_22 = -one(q)

    σ_z = (σ_z_11 ⊕ σ_z_22)/2

    H_kinetic = op_px^2/(2*m) + op_py^2/(2*m)

    σ_px = (op_px ⊕ (-op_px))/2

    σ_px_I = (one(p) ⊕ one(p))/4

    #spinor_H_kinetic = (H_kinetic ⊕ H_kinetic) -1im*(1/m_dimless)*k*σ_pz + (1/(2*m_dimless))*k^2*σ_pz_I
    spinor_H_kinetic = (H_kinetic ⊕ H_kinetic) - (1/m)*k*σ_px + (1/(2*m))*(k^2)*σ_px_I
    #spinor_H_kinetic = (H_kinetic ⊕ H_kinetic) - (1/m_dimless)*k_dimless*σ_pz

    spinor_H_kinetic_FFT = LazyProduct(spinor_Tqp,spinor_H_kinetic,spinor_Tpq)

    H1 = (Ω/2)*σ_x + (δ/2)*σ_z

    return LazySum(spinor_H_kinetic_FFT,V,H1,Hψ)

end

H_total = makeHamiltonian(b_position_x,b_position_y,b_position,b_momentum_x,b_momentum_y,b_momentum,m_dimless,k_dimless,Ω,δ,V_trap,Hψ)

function H_im(t, ψ) # Update state-dependent term in H
    Hψ.data.nzval[1:steps] .= g_uu/(dx*dy) * abs2.(getblock(ψ,1).data) + g_ud/(dx*dy) * abs2.(getblock(ψ,2).data)
    Hψ.data.nzval[steps+1:end] .= g_dd/(dx*dy) * abs2.(getblock(ψ,2).data) + g_du/(dx*dy) * abs2.(getblock(ψ,1).data)
    return H_total
end

spinor_ψ0 = initialState(b_position_x,b_position_y,x0,y0,p0_x,p0_y,σ,separation)

# renormalization callback
norm_func(u, t, integrator) = normalize!(u)
ncb = FunctionCallingCallback(norm_func; func_everystep = true, func_start = true)

tout, ψt = timeevolution.schroedinger_dynamic(T, spinor_ψ0, H_im, callback=ncb, alg_hints=[:stiff], maxiters=maxiters_int)

density_up = [Array(transpose(reshape((abs2.(getblock(ψ,1).data)), (Npoints_x, Npoints_y)))) for ψ = ψt]
density_down = [Array(transpose(reshape((abs2.(getblock(ψ,2).data)), (Npoints_x, Npoints_y)))) for ψ = ψt]

xsample, ysample = samplepoints(b_position_x), samplepoints(b_position_y)

jldsave(string("diamond.jld2");xsample,ysample,tout,ψt,Npoints_x,Npoints_y)

plotDiamond("diamond.png",V_plot,5)
