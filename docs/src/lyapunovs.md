The Finite Time Lyapunov Spectrum (FTLS) for a 2D billiard system consists of a
set of 4 numbers $\lambda_i \, , \{ i = 1, ...,4 \}$ that characterize how fast
the separation of initially close initial conditions grows.

It can be shown theoretically that two of these exponents must be zero
($\lambda_2$ =$\lambda_3$ = 0) and the other two are paired in such a way that
they sum up to zero, i.e. $\lambda_1 =  -\lambda_4$).

The function provided to calculate the FTLS is
```@docs
lyapunovspectrum
```

Here its basic use is illustrated
```@example lyaps
using DynamicalBilliards

radius = 1.0
l = 2.0

bd = Billiard(billiard_polygon(6, l; setting = "periodic")..., Disk([0., 0.], radius))

par = randominside(bd)
t = 1000.0

exps = lyapunovspectrum(par, bd, t)
```

In the following example we compute the change of $\lambda_1\$ versus the
distance between the disks in a hexagonal periodic billiard.

```@example lyaps
using DynamicalBilliards
using PyPlot

t = 5000.0
radius = 1.0

spaces = 2.0:0.1:4.4 #Distances between adjacent disks
lyap_time = zero(spaces) #Array where the exponents will be stored

for (i, space) in enumerate(spaces)
    bd = billiard_polygon(6, space/(sqrt(3)); setting = "periodic")
    disc = Disk([0., 0.], radius)
    billiard = Billiard(bd.obstacles..., disc)
    p = randominside(billiard)
    lyap_time[i] = lyapunovspectrum(p, billiard, t)[1]
end
figure()
plot(spaces, lyap_time, "*-")
xlabel("\$w\$"); ylabel("\$\\lambda_1\$")
savefig("lyapos.svg"); nothing # hide
```
![](lyapos.svg)

The plot of the maximum exponent can be compared with the results reported by
[Gaspard et. al](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.51.5332)
(see figure 7), showing that using just `t = 5000.0` is already enough of a
statistical averaging.

## Perturbation Growth
To be able to inspect the dynamics of perturbation growth in more detail, we also provide the following function:
```@docs
perturbationgrowth
```
---
For example, lets plot the evolution of the perturbation growth using different colors for collisions with walls and disks in the Sinai billiard:
```@example lyaps
using DynamicalBilliards, PyPlot, LinearAlgebra
bd = billiard_sinai()

ts, Rs, is = perturbationgrowth(Particle(0.1, 0.1, 0.1), bd, 10.0)
Δ = perturbationevolution(Rs)

figure()
plot(ts, log.(norm.(Δ)), "k-", lw = 0.5)
scatter(ts, log.(norm.(Δ)), c = [j == 1 ? "C0" : "C1" for j in is])
xlabel("\$t\$"); ylabel("\$\\log(||\\Delta ||)\$")
savefig("pertg.svg"); nothing # hide
```
![](pertg.svg)
