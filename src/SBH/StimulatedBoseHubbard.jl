module StimulatedBoseHubbard

export LatticeBasis
export initial, rand
export latticeSolve
export averages

using QuantumOptics
using QuantumOpticsBase
using SparseArrays
using OrdinaryDiffEq, DiffEqCallbacks
using LinearAlgebra

using Unzip

import Base: ==, ones, rand

import LinearAlgebra: normalize!

"""
    ActionBasis(P)

Basis for a Action space where `P` specifies an action.
"""
struct ActionBasis{T} <: Basis
    shape::Vector{T}
    N::T
    function ActionBasis(N::T) where T
        new{T}([N+1], N)
    end
end

==(b1::ActionBasis, b2::ActionBasis) = (b1.N==b2.N)

function uniform(b::ActionBasis)
    normalize!(Ket(b,ones(b.N)))
end

function initial(b::ActionBasis,n)
    data = ComplexF64.(zeros(b.N+1))
    data[n+1] = ComplexF64(exp(im*rand()*2*π))
    Ket(b,data)
end

"""
    action([T=ComplexF64,] b::ActionBasis)

Action operator for the specified Action space with optional data type `T`.
"""
function action(::Type{C}, b::ActionBasis) where C
    T = real(C)
    diag = @. C(T(0:b.N))
    data = spdiagm(0 => diag)
    SparseOperator(b, data)
end

action(b::ActionBasis) = action(ComplexF64, b)

"""
    sq_action([T=ComplexF64,] b::ActionBasis)

square root of the Action operator for the specified
    Action space with optional data type `T`.
"""
function sq_action(::Type{C}, b::ActionBasis) where C
    T = real(C)
    diag = @. C(sqrt(T(0:b.N)))
    data = spdiagm(0 => diag)
    SparseOperator(b, data)
end

sq_action(b::ActionBasis) = sq_action(ComplexF64, b)


"""
    destroy([T=ComplexF64,] b::ActionBasis)

Annihilation operator for the specified Action space with optional data type `T`.
"""
function destroy(::Type{C}, b::ActionBasis) where C
    T = real(C)
    diag = @. C(sqrt(T(1:b.N)))
    data = spdiagm(1 => diag)
    SparseOperator(b, data)
end

destroy(b::ActionBasis) = destroy(ComplexF64, b)


"""
    create([T=ComplexF64,] b::ActionBasis)

Creation operator for the specified Action space with optional data type `T`.
"""
function create(::Type{C}, b::ActionBasis) where C
    T = real(C)
    diag = @. C(sqrt(T(1:b.N)))
    data = spdiagm(-1 => diag)
    SparseOperator(b, data)
end

create(b::ActionBasis) = create(ComplexF64, b)

"""

Create the identity operator for Action space
"""
function ones(::Type{C}, b::ActionBasis) where C
    identityoperator(C,b)
end

ones(b::ActionBasis) = ones(ComplexF64, b)

"""

Create the zero operator for Action space
"""
function zero(::Type{C}, b::ActionBasis) where C
    l = b.N+1
    data = spzeros(l,l)
    SparseOperator(b,data)
end

zero(b::ActionBasis) = zero(ComplexF64, b)

"""
    AngleBasis(Q)

Basis for a Angle space wwhere 'Q' specifies an angle.
"""
struct AngleBasis{T} <: Basis
    shape::Vector{T}
    N::T
    function AngleBasis(N::T) where T
        new{T}([N], N)
    end
end

==(b1::AngleBasis, b2::AngleBasis) = (b1.N==b2.N)

function initial(b::AngleBasis)
    p = rand(1:b.N)
    q = zeros(ComplexF64,b.N)
    q[p] = exp(im*rand()*2*π)
    Ket(b,q)
end

"""
    phase([T=ComplexF64,] b::AngleBasis)

Phase operator for the specified Angle space with optional data type `T`.
"""
function phase(::Type{C}, b::AngleBasis) where C
    T = real(C)
    diag = @. C(T((0:(b.N-1)).*(2*π./b.N)))
    data = spdiagm(0 => diag)
    SparseOperator(b, data)
end

phase(b::AngleBasis) = phase(ComplexF64, b)

"""
    shift([T=ComplexF64,] b::AngleBasis)

Operator that shifts phase with optional type `T`.
"""
function shift(::Type{C}, b::AngleBasis, ϕ::Float64) where C
    T = real(C)
    if ϕ >= 2*π*(1/b.N)
        p = -1
    elseif ϕ <= -2*π*(1/b.N)
        p = 1
    else
        p = 0
    end
    diag = C.(ones(b.N-abs(p)))
    data = spdiagm(p => diag)
    if p == -1
        data[1,b.N] = 1
    elseif p == 1
        data[b.N,1] = 1
    end
    SparseOperator(b, data)
end

shift(b::AngleBasis, ϕ) = shift(ComplexF64, b, ϕ)

"""
    ones(b)

Create the identity operator for Angle space
"""
function ones(::Type{C}, b::AngleBasis) where C
    identityoperator(C,b)
end

ones(b::AngleBasis) = ones(ComplexF64, b)

"""

Create the zero operator for Angle space
"""
function zero(::Type{C}, b::AngleBasis) where C
    data = spzeros(l,l)
    SparseOperator(b,data)
end

zero(b::AngleBasis) = zero(ComplexF64, b)

"""
    ActionAngleBasis(NP,NQ)

Basis for ActionAngle space where `NP` is the number of actions and
    `NQ` is granularity of phase
"""
struct ActionAngleBasis{T} <: Basis
    shape::Vector{T}
    NP::T
    NQ::T
    function ActionAngleBasis(NP::T,NQ::T) where T
        new{T}([(NP+1)+NQ], NP, NQ)
    end
end

function initial(b::ActionAngleBasis,n)
    b1 = initial(ActionBasis(b.NP),n)
    b2 = initial(AngleBasis(b.NQ))

    Ket(b,[b1.data;b2.data])
end

function makeActionAngleOperator(::Type{C}, b::ActionAngleBasis, f, g) where C
    o1 = f(ActionBasis(b.NP))
    o2 = g(AngleBasis(b.NQ))

    SparseOperator(b,blockdiag(o1.data,o2.data))
end

function makeActionAngleOperator(::Type{C}, b::ActionAngleBasis, f, g, ϕ) where C
    o1 = f(ActionBasis(b.NP))
    o2 = g(AngleBasis(b.NQ),ϕ)

    SparseOperator(b,blockdiag(o1.data,o2.data))
end

function action(::Type{C}, b::ActionAngleBasis) where C
    makeActionAngleOperator(C,b,action,ones)
end

action(b::ActionAngleBasis) = action(ComplexF64, b)

function sq_action(::Type{C}, b::ActionAngleBasis) where C
    makeActionAngleOperator(C,b,sq_action,ones)
end

sq_action(b::ActionAngleBasis) = sq_action(ComplexF64, b)

function destroy(::Type{C}, b::ActionAngleBasis) where C
    makeActionAngleOperator(C,b,destroy,ones)
end

destroy(b::ActionAngleBasis) = destroy(ComplexF64, b)

function create(::Type{C}, b::ActionAngleBasis) where C
    makeActionAngleOperator(C,b,create,ones)
end

create(b::ActionAngleBasis) = create(ComplexF64, b)

function phase(::Type{C}, b::ActionAngleBasis) where C
    makeActionAngleOperator(C,b,ones,phase)
end

phase(b::ActionAngleBasis) = phase(ComplexF64, b)

function shift(::Type{C}, b::ActionAngleBasis, ϕ) where C
    makeActionAngleOperator(C,b,ones,shift,ϕ)
end

shift(b::ActionAngleBasis, ϕ) = shift(ComplexF64, b, ϕ)

function ones(::Type{C}, b::ActionAngleBasis) where C
    identityoperator(C,b)
end

ones(b::ActionAngleBasis) = ones(ComplexF64, b)

function zero(::Type{C}, b::ActionAngleBasis) where C
    makeActionAngleOperator(C,b,zero,zero)
end

zero(b::ActionAngleBasis) = zero(ComplexF64, b)

function aadiag(::Type{C}, b::ActionAngleBasis, i, k) where C
    data = spdiagm(i => k.data)
    SparseOperator(b, data)
end

aadiag(b::ActionAngleBasis, i, Ψ) = aadiag(ComplexF64, b, i, Ψ)

struct LatticeBasis{T} <: Basis
    shape::Vector{T}
    N::T
    NP::T
    NQ::T
    function LatticeBasis(N::T,NP::T,NQ::T) where T
        new{T}([N*((NP+1)+NQ)], N, NP, NQ)
    end
end

==(b1::LatticeBasis, b2::LatticeBasis) = (b1.N == b2.N && b1.NP==b2.NP && b2.NQ == b2.NQ)

function initial(l::LatticeBasis,n::Integer)
    b = ActionAngleBasis(l.NP,l.NQ)

    ket = initial(b,n)

    d = ket.data

    for i = 2:l.N
        ket = initial(b,n)

        d = vcat(d,ket.data)
    end

    return Ket(l,d)
end

function rand(l::LatticeBasis,n::Integer)
    b = ActionAngleBasis(l.NP,l.NQ)

    m = rand(0:n,l.NP+1)

    m = (m ./ sum(m)) * n

    m = Int64.(round.(m))

    ket = initial(b,m[1])

    d = ket.data

    for i = 2:l.N
        ket = initial(b,m[i])

        d = vcat(d,ket.data)
    end

    return Ket(l,d)
end

function averages(l::LatticeBasis,Ψ)
    data = Ψ.data

    b1 = ActionBasis(l.NP)
    b2 = AngleBasis(l.NQ)

    d1 = l.NP+1
    d2 = l.NQ
    d = d1 + d2

    as = []

    for i = 1:l.N
        aa = data[((i-1)*d+1):i*d]
        a1 = aa[1:d1]
        a2 = aa[(d1+1):d]

        n = expect(action(b1),Ket(b1,a1))
        t = expect(phase(b2),Ket(b2,a2))

        n_v = variance(action(b1),Ket(b1,a1))
        t_v = variance(phase(b2),Ket(b2,a2))

        push!(as,(real(n),real(t),real(n_v),real(t_v)))
    end

    return as
end

function normalize!(data,l::LatticeBasis)

    b1 = ActionBasis(l.NP)
    b2 = AngleBasis(l.NQ)

    d1 = l.NP+1
    d2 = l.NQ
    d = d1 + d2

    aa = data[1:d]
    a1 = aa[1:d1]
    a2 = aa[(d1+1):d]

    n = Ket(b1,a1)
    t = Ket(b2,a2)

    normalize!(n)
    normalize!(t)

    dd = vcat(n.data,t.data)

    for i = 2:l.N
        aa = data[((i-1)*d+1):i*d]
        a1 = aa[1:d1]
        a2 = aa[(d1+1):d]

        n1 = Ket(b1,a1)
        t1 = Ket(b2,a2)

        normalize!(n1)
        normalize!(t1)

        dd = vcat(dd,vcat(n1.data,t1.data))
    end

    for i = 1:(l.N*d)
        data[i] = dd[i]
    end

end

function ones(::Type{C},l::LatticeBasis) where C
    identityoperator(C,l)
end

ones(l::LatticeBasis) = ones(ComplexF64,l)

function makeLatticeOperator(::Type{C}, l::LatticeBasis, f) where C
    b = ActionAngleBasis(l.NP,l.NQ)
    op = f(b)

    d = op.data

    for i = 2:l.N
        op1 = f(b)

        d = blockdiag(d,op1.data)
    end
    SparseOperator(l,d)
end

makeLatticeOperator(l::LatticeBasis, f) = makeLatticeOperator(ComplexF64, l, f)

function makeLatticeOperator(::Type{C}, l::LatticeBasis, f, a) where C
    b = ActionAngleBasis(l.NP,l.NQ)
    op = f(b,a)

    d = op.data

    for i = 2:l.N
        op1 = f(b,a)

        d = blockdiag(d,op1.data)
    end
    SparseOperator(l,d)
end

makeLatticeOperator(l::LatticeBasis, f, a) = makeLatticeOperator(ComplexF64, l, f, a)

function makeLatticeOperatorIndex(::Type{C}, l::LatticeBasis, i, f) where C
    b = ActionAngleBasis(l.NP,l.NQ)
    op = f(b)

    if i == 1
        d = op.data
    else
        d = ones(b).data
    end

    for k = 2:l.N
        if i == k
            d = blockdiag(d,op.data)
        else
            d = blockdiag(d,ones(b).data)
        end
    end
    SparseOperator(l,d)
end

makeLatticeOperatorIndex(l::LatticeBasis, i, f) = makeLatticeOperatorIndex(ComplexF64, l, i, f)

function makeLatticeOperatorIndex(::Type{C}, l::LatticeBasis, i, f, a) where C
    b = ActionAngleBasis(l.NP,l.NQ)
    op = f(b,a)

    if i == 1
        d = op.data
    else
        d = ones(b).data
    end

    for k = 2:l.N
        if i == k
            d = blockdiag(d,op.data)
        else
            d = blockdiag(d,ones(b).data)
        end
    end
    SparseOperator(l,d)
end

makeLatticeOperatorIndex(l::LatticeBasis, i, f, a) = makeLatticeOperatorIndex(ComplexF64, l, i, f, a)

function action(::Type{C}, l::LatticeBasis) where C
    makeLatticeOperator(C,l,action)
end

action(l::LatticeBasis) = action(ComplexF64, l)

function sq_action(::Type{C}, l::LatticeBasis) where C
    makeLatticeOperator(C,l,sq_action)
end

sq_action(l::LatticeBasis) = sq_action(ComplexF64, l)

function destroy(::Type{C}, l::LatticeBasis) where C
    makeLatticeOperator(C,l,destroy)
end

destroy(l::LatticeBasis) = destroy(ComplexF64, l)

function create(::Type{C}, l::LatticeBasis) where C
    makeLatticeOperator(C,l,create)
end

create(l::LatticeBasis) = create(ComplexF64, l)

function phase(::Type{C}, l::LatticeBasis) where C
    makeLatticeOperator(C,l,phase)
end

phase(l::LatticeBasis) = phase(ComplexF64, l)

function shift(::Type{C}, l::LatticeBasis, ϕ) where C
    makeLatticeOperator(C,l,shift,ϕ)
end

shift(l::LatticeBasis, ϕ) = shift(ComplexF64, l, ϕ)

function extractBlock(m,l::LatticeBasis,i,j)
    I, J, V = findnz(m)
    d = l.NP+1+l.NQ

    II = findall(x -> x > (i-1)*d && x <= i*d,I)
    JJ = findall(x -> x > (j-1)*d && x <= j*d,J)

    K = intersect(II,JJ)

    III = I[K] .- (i-1)*d
    JJJ = J[K] .- (j-1)*d

    VVV = V[K]

    n = sparse(III,JJJ,VVV,d,d)

    return n
end

function addBlockToHamiltonian(H,l::LatticeBasis,i,j,op)
    h = H.data
    o = op.data

    p = l.NP+1+l.NQ

    (hI,hJ,hV) = findnz(h)
    (oI,oJ,oV) = findnz(o)

    updated = false

    for m in 1:length(oI)
        nn = findall(x -> x == oI[m]+p*(i-1),hI)
        updated = false
        for n in nn
            hn = hJ[n]
            if hn == oJ[m]+p*(j-1)
                H.data[hI[n],hn] += oV[m]
                updated = true
            end
        end
        if updated == false
            H.data[oI[m]+p*(i-1),oJ[m]+p*(j-1)] = oV[m]
        end
    end

    return H
end

function replaceBlockInHamiltonian(H,l::LatticeBasis,i,j,op)
    h = H.data
    o = op.data

    p = l.NP+1+l.NQ

    (hI,hJ,hV) = findnz(h)

    I = findall(x -> x > (i-1)*p && x <= i*p,hI)
    J = findall(x -> x > (j-1)*p && x <= j*p,hJ)

    K = intersect(I,J)

    for k in K
        h[hI[k],hJ[k]] = 0
    end

    (oI,oJ,oV) = findnz(o)

    for k in 1:length(oI)
        H.data[oI[k]+(i-1)*p,oJ[k]+(j-1)*p] = oV[k]
    end

    dropzeros!(H.data)

    return H
end

function makeHamiltonian(::Type{C},l::LatticeBasis,K,α,β,ϵ,U,J,φ,Ψ) where C
    basis = ActionAngleBasis(l.NP,l.NQ)

    as = averages(l,Ψ)

    H_cp = -ϵ*action(C,l)
    H_ons = (U/2)*action(C,l)*(action(C,l)-ones(l))

    H_int = ones(l)

    for i = 1:l.N
        for j = 1:l.N
            if K[i,j] != 0

                ci = makeLatticeOperatorIndex(l,i,create)
                dj = makeLatticeOperatorIndex(l,j,destroy)

                dj = addBlockToHamiltonian(dj,l,i,j,destroy(C,basis))

                ni = makeLatticeOperatorIndex(l,i,action)

                ni = addBlockToHamiltonian(ni,l,i,j,-action(C,basis))

                nj = makeLatticeOperatorIndex(l,j,action)

                nj = addBlockToHamiltonian(nj,l,j,i,-action(C,basis))

                ϕ = as[i][2] - as[j][2]

                s = makeLatticeOperatorIndex(l,i,shift,ϕ/2)

                H_int += -K[i,j]*J*ci*dj*((1-α)*ones(l)+α*(ni-nj+ones(l)))*((1-β)*ones(l)+β*s*cos(ϕ+φ))
            end
        end
    end
    return H_cp + H_ons + H_int
end

makeHamiltonian(l::LatticeBasis,c,α,β,ϵ,U,J,φ,Ψ) = makeHamiltonian(ComplexF64,l,c,α,β,ϵ,U,J,φ,Ψ)

function latticeSolve(p,Ψ0)
    T = p["T"]
    b = LatticeBasis(p["nodes"],p["particles"],p["angles"])
    c = p["edges"]
    α = p["α"]
    β = p["β"]
    ϵ = p["ϵ"]
    U = p["U"]
    J = p["J"]
    φ = p["φ"]

    H_sbh = (t,Ψ) -> makeHamiltonian(b,c,α,β,ϵ,U,J,φ,Ψ)

    norm_func = (u, t, integrator) -> normalize!(u,b)

    ncb = FunctionCallingCallback(norm_func; func_everystep = true, func_start = true)

    maxiters_int = 1e8

    @time tout, Ψt = timeevolution.schroedinger_dynamic(T, Ψ0, H_sbh, callback=ncb, alg_hints=[:stiff], maxiters=maxiters_int)

    return tout, Ψt
end

end
