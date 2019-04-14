using DynamicalBilliards, PyPlot
using DynamicalBilliards: normalvec, cossin

# superellipse curve
function superellipse(a, b, n)
    ρ(φ) = 10.0/sqrt(sqrt((5*cos(φ))^4 + (2*sin(φ))^4))
end

ρ(φ) = 10.0/sqrt(sqrt((5*cos(φ))^4 + (2*sin(φ))^4))
w(x) = -2.5*(64*sin(x)^3*cos(x) - 2500*sin(x)*cos(x)^3) /
        ((5*cos(x))^4 + (2*sin(x))^4)^1.25

l = PolarCurve([0.0, 0.0], ρ, w, "Superellipse-4")

figure()

plot(l)

for φ in 0:0.1:2π
    ρ = l.ρ(φ)
    co, si = cossin(φ)
    n = normalvec(l, ρ .* cossin(φ))
    plot([ρ*co, ρ*co + n[1]], [ρ*si, ρ*si + n[2]], color = "r")
end
