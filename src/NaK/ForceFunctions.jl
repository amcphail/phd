
module ForceFunctions

export F_spring
export F_pumping, F_forcing, F_damping
export F_lattice

function F_spring(z,k=1)
    return - k * z
end

# Used maxima to integrate
# z: sheath middle
# c: core middle
# z: core radius
#
# e1 : ((z-c)/((z-c)^2+a^2)^(3/2));
#
# b: half height
# d: sheath middle
#
# integrate(integrate(e1,c,-b,b),z,d-b,d+b);
#
function F_pumping(z,a=1,b=1,k=1,rho_c=1,rho_s=1)
    integral = - asinh((z+2*b)/a) - asinh((z-2*b)/a) + 2*asinh(z/a)
    return - 2*pi*a*k*rho_c*rho_s * integral
end

function F_forcing(z,armed,ran,k=1.,p=0.01)
    if armed < eps(0.0)
        if z < 0
            if ran < p
                return - k
            end
        end
    end
    return 0
end

function F_damping(v,k=1.)
    return - k * v
end

function F_lattice(z1,z2,a=1.,k=10)
    θ = atan(z1-z2,a)
    r = sqrt(a^2+(z2-z1)^2)

    return -k*2*sin(θ) / r
end

end
