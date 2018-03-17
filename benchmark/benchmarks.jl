using BenchmarkTools, DynamicalBilliards, IterTools


const SUITE = BenchmarkGroup(["DynamicalBilliards"])
bt = billiard_mushroom()
bt2 = billiard_sinai(;setting="periodic")
particles = [Particle(0.05, 0.05, -0.1), MagneticParticle(0.05,0.05,-0.1,1.0)]
ptypes = ["straight", "magnetic"]
colf = (collisiontime, next_collision)
name = (f) -> split(string(f), '.')[end]

for f in colf
    SUITE[name(f)] = BenchmarkGroup(["propagation", "collision"])
    for ptype ∈ ptypes
        SUITE[name(f)][ptype] = BenchmarkGroup(["propagation", "collision", ptype])
    end
end

for (f, p) in zip(["straight", "magnetic"], particles)
    for o in chain(bt, bt2)
        SUITE["collisiontime"][f][o.name] =
            @benchmarkable collisiontime($p, $o)
    end
end

for (f, p) in zip(["straight", "magnetic"], particles)
    for (bname, bil) in zip(["mushroom", "psinai"], (bt, bt2))
        SUITE["next_collision"][f][bname] =
            @benchmarkable next_collision($p, $bil)
    end
end
