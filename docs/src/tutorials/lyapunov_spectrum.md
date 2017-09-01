The Finite Time Lyapunov Spectrum (FTLS) for the billiard dynamics consists of a set of 4 numbers $\lambda_i \, , \{ i = 1, ...,4 \}$ that characterizes the chaoticity of the billiard. They depend on the initial condition and the total time of integration.

It can be shown theoretically that two of these exponents must be zero ($\lambda_2 = \lambda_3 = 0$) and the other two are paired in such a way that they sum up to zero, i.e. $\lambda_1 =  -\lambda_4$).

The function provided to calculate the FTLS is
```
 lyapunovspectrum(p::Particle, bt::Vector{Obstacle}, t::Float64)
```
and it returns an array with the 4 lyapunov exponents.

Here its basic use is illustrated
```
radius = 1.0
l = 1.5

bt = billiard_polygon(6, 1.5; setting = "periodic")
disc = Disk([0., 0.], radius)
push!(bt, disc)

p = randominside(bt)
t = 1000.0

exps = lyapunovspectrum(p, bt, t)
```

The following code is for a family of polygonal billiards (hexagonal unit cell) parameterized by the space between the disks.

```
using DynamicalBilliards
radius = 1.0
space = 3.0
sides = 6
t = 5000.0

sp = 2.1:0.1:4.4
lyap_time = zeros(sp)

disc = Disk([0., 0.], radius)

i=1
for space in sp
	bt = billiard_polygon(sides, space/(2*cos(pi/sides)); setting = "periodic")
	push!(bt, disc)
    p = randominside(bt)
    exps = lyapunovspectrum(p, bt, t)
    lyap_time[i] = exps[1]
    i+=1
end

```

The plot of the maximum exponent is displayed below and it can be compared with the results reported by [Gaspard et. al](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.51.5332)(figure 7. also displayed here) for the average over an ensemble.

![Lyapunov Exponent](lyap.png)
![gaspard](gaspard.png)


