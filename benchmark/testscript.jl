using Revise; using DynamicalBilliards, BenchmarkTools

bt = billiard_sinai()
p = MagneticParticle(0.1,0.1, π/8, 0.5)

@btime next_collision($p, $bt)

  197.933 ns (0 allocations: 0 bytes)
(0.3561834976110847, 1)
