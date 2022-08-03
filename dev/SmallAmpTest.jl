using AKNS, ApproxFun, ApproxFunRational, FastGaussQuadrature, Plots, SpecialFunctions

function remove_deltas!(s::AKNS.SumFun)
    i = 0
    while i < length(s.funs)
        i += 1
        if typeof(s.funs[i].space) <: AKNS.DiracSpace
            deleteat!(s.funs,i)
            i -= 1
        end
    end
    return s
end

## Explicit reflection test for focusing case
A = .01; γ = 1.0;
q = x -> -1im*A*sech(x)*exp(-1im*γ*A*log(cosh(x)))
r = x -> -conj(q(x))
w = z -> -1im*z - A*γ*1im/2 + 0.5
T = sqrt(γ^2/4 -1 |> Complex)
wp = -1im*A*(T + γ/2)
wm = 1im*A*(T-γ/2)
asech = z -> (gamma(w(z))*gamma(w(z)-wm-wp))/(gamma(w(z)-wp)*gamma(w(z)-wm))
bsech = z -> 1im*1/A*2^(-1im*A*γ)*(gamma(w(z))*gamma(1-w(z)+wp+wm))/(gamma(wp)*gamma(wm))
ρsech = z -> bsech(z)/asech(z)
logasech = z -> loggamma(w(z)) + loggamma(w(z)-wm-wp) - loggamma(w(z)-wp) - loggamma(w(z)-wm)

out = AKNSscattering(q,r,200,400)
function ρ_notfun(x,out)
    s = out(x)
    return s[1,2]/s[1,1]
    #return s[2,1]/s[1,1]  # = ρsech
end

F = FourierTransform(-1.0)
## Fourier transform of potential
qhat = F*Fun(zai(q),OscLaurent(0.0,1.0),200)

x = -5:0.011:5
# FT of potential on grid
y1 = map(x -> qhat(x),x)
# True reflection coefficient
y2 = map(x -> ρsech(x/2),x)
# Numerically computed reflection coefficient
y3 = map(x -> ρ_notfun(x/2,out),x)

ρsech(.3)
ρ_notfun(.3,out)

plot(x,abs.(y1/A))
plot!(x,abs.(y2/A))
plot!(x,abs.(y3/A))

plot(x,real(y1/A))
plot!(x,real(y2/A))
plot!(x,real(y3/A))

plot(x,imag(y1/A))
plot!(x,imag(y2/A))
plot!(x,imag(y3/A))

## Reflection test for defocusing case
A = .01; γ = 1.0;
q = x -> -1im*A*x*sech(x)*exp(-1im*γ*A*log(cosh(x)))
r = x -> conj(q(x))

out = AKNSscattering(q,r,200,400)
function ρ_notfun(x,out)
    s = out(x)
    return s[1,2]/s[1,1]
    #return s[1,2]/s[2,2]
end

F = FourierTransform(-1.0)
## Fourier transform of potential
qhat = F*Fun(zai(q),OscLaurent(0.0,1.0),200)

x = -5:0.011:5
# FT of potential on grid
y1 = map(x -> qhat(x),x)
# Numerically computed reflection coefficient
y2 = map(x -> ρ_notfun(x/2,out),x)

## Compare FT of poential with reflection
plot(x,abs.(y1/A))
plot!(x,abs.(y2/A))

plot(x,real(y1/A))
plot!(x,real(y2/A))

plot(x,imag(y1/A))
plot!(x,imag(y2/A))

F = FourierTransform(1.0) # negative in fourier exponent

#### THIS IS THE KERNEL FUNCTION ####
ρ = 1/(2*pi)*(F*Fun(zai(x -> ρ_notfun(x/2,out)),OscLaurent(0.0,1.0),200)) |> AKNS.SumFun |> remove_deltas!;

x = -5:0.011:5
# Potential on grid
y1 = map(x -> q(x),x)
# FT of computed reflection coefficient
y2 = map(ρ,x)

## Compare poential with FT of reflection
plot(x,abs.(y1/A))
plot!(x,abs.(y2/A))

plot(x,real(y1/A))
plot!(x,real(y2/A))

plot(x,imag(y1/A))
plot!(x,imag(y2/A))