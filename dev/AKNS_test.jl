using ApproxFunRational, AKNS, SpecialFunctions, ApproxFun
using AbstractIterativeSolvers

A = 0.65;
γ = 2.0;
w = z -> -1im * z - A * γ * 1im / 2 + 0.5
T = sqrt(γ^2 / 4 - 1 |> Complex)
wp = -1im * A * (T + γ / 2)
wm = 1im * A * (T - γ / 2)

q = x -> -1im * A * sech(x) * exp(-1im * γ * A * log(cosh(x)))
r = x -> -conj(q(x))

logasech = z -> loggamma(w(z)) + loggamma(w(z) - wm - wp) - loggamma(w(z) - wp) - loggamma(w(z) - wm)
logbsech = z -> loggamma(w(z)) + loggamma(1 - w(z) + wp + wm) - loggamma(wp) - loggamma(wm)
ρsechexplog = z -> 1im * 1 / A * 2^(-1im * γ * A) * exp(logbsech(z) - logasech(z))
bsechexplog = z -> 1im * 1 / A * 2^(-1im * γ * A) * exp(logbsech(z))
asechexplog = z -> exp(logasech(z))


AD = AKNSdet(z -> ρsechexplog(z / 2) / (2 * pi), z -> -conj(ρsechexplog(z / 2)) / (2 * pi), 60, 300)

x = 0.1

AD.f1(20)

function potential(x)
    M = abs(x) + 10
    xs = M * (AD.quad[1] .+ 1) / 2
    ω = sqrt.(AD.quad[2] * M / 2)
    W = Diagonal(ω)
    f = (X, Y) -> AD.f1(x + X + Y + 5e-16)
    Ax = W * f.(xs, xs') * W
    f = (X, Y) -> AD.f2(x + X + Y + 5e-16)
    Bx = W * f.(xs, xs') * W
    Avec = W * map(z -> AD.f1(x + z + 5e-16), xs)
    AD.f1(x) - Avec' * ((I + Bx * Ax) \ (Bx * Avec))
end

potential(0.5)

x = 0:0.1:8
y = map(potential, x)

using Plots

plot(x, real(y), legend=false)
plot!(x, imag(y))
plot!(x, real(map(q, x)))
plot!(x, imag(map(q, x)))
