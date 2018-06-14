using Revise; using DynamicalBilliards, BenchmarkTools

bt = billiard_sinai()
p = MagneticParticle(0.1,0.1, π/8, 0.5)

@btime next_collision($p, $bt)

@btime bounce!($p, $bt)



# before sincos():
julia> @btime bounce!($p, $bt)
  278.529 ns (0 allocations: 0 bytes)
(1, 0.44691013238163385, [0.707091, 0.640047], [-0.353179, 0.935556])

# after sincos():
julia> @btime bounce!($p, $bt)
  57.476 ns (0 allocations: 0 bytes)
(0, Inf, [0.415228, 1.10629], [0.957396, -0.288778])
