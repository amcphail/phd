height = round(30e-9/E_k) # nK

function wells(x, y, d)
    if ((x > -4 && x < 0) || ((x > d) && (x < (4 + d)))) && (y > -2 && y < 2)
        return 0
    else
        return height
    end
end

distance = 1

potential(x,y) = wells(x,y,distance)

V = potentialoperator(comp_q, potential)

Hg = diagonaloperator(b_position_x, Ket(b_position_x).data)⊗diagonaloperator(b_position_y, Ket(b_position_y).data) # ∝ |ψ|^2

H = -1im*LazySum(H_kin_x_FFT,H_kin_y_FFT,V,Hg)

dx = (xmax-xmin)/Npoints_x
dy = (ymax-ymin)/Npoints_y

function Hgp(t, ψ) # Update state-dependent term in H
    normalize!(ψ)
    H.operators[4].data.nzval .= g_2D*abs2.(ψ.data)/(dx*dy)
    return H
end



T = collect(0.0:1:10)*0.5

tout, ψt = timeevolution.schroedinger_dynamic(T, ψ, Hgp)
