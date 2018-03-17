using BenchmarkTools, DynamicalBilliards, IterTools


const SUITE = BenchmarkGroup(["DynamicalBilliards"])
bt = billiard_mushroom()
bt2 = billiard_sinai()
particles = [Particle(0.05, 0.05, -0.1), MagneticParticle(0.05,0.05,-0.1,1.0)]
ptypes = ["straight", "magnetic"]


SUITE["coltimes"] = BenchmarkGroup(["propagation", "collision"])
for ptype ∈ ptypes
    SUITE["coltimes"][ptype] = BenchmarkGroup(["propagation", "collision", ptype])
    for obstype ∈ ["obstacles", "biltable"]
        SUITE["coltimes"][ptype][obstype] =
            BenchmarkGroup(["propagation", "collision", ptype])
    end
end

for (f, p) in zip(["straight", "magnetic"], particles)
    for o in chain(bt, bt2)
        SUITE["coltimes"][f]["obstacles"][o.name] =
            @benchmarkable collisiontime($p, $o)
    end
end

# SUITE["trigonometry"]["hyperbolic"] = BenchmarkGroup()
# for f in (sin, cos, tan)
#     for x in (0.0, pi)
#         SUITE["trigonometry"]["hyperbolic"][string(f), x] = @benchmarkable ($f)($x)
#     end
# end

# showall(SUITE)
