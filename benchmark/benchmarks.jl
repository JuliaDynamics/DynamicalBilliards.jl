using BenchmarkTools, DynamicalBilliards, IterTools


const SUITE = BenchmarkGroup(["DynamicalBilliards"])
bt = billiard_mushroom()
bt2 = billiard_sinai(;setting="periodic")

################################################################################
## COLLISION & PROPAGATION
################################################################################

particles = [Particle(0.05, 0.05, -0.1), MagneticParticle(0.05,0.05,-0.1,1.0)]
obstacles = [bt[1], bt[6], bt2[1], bt2[5]] #distinct obstacles for resolvecollision! tests
proptime = 4.2
ptypes = ["straight", "magnetic"]
colf = (collision,
        next_collision,
        bounce!,
        resolvecollision!,
        propagate!,
        lyapunovspectrum
        )
name = (f) -> split(string(f), '.')[end]

for f in colf
    SUITE[name(f)] = BenchmarkGroup(["propagation", "collision"])
    for ptype ∈ ptypes
        SUITE[name(f)][ptype] = BenchmarkGroup(["propagation", "collision", ptype])
    end
end

for (f, p) in zip(["straight", "magnetic"], particles)
    for o in chain(bt, bt2)
        SUITE["collision"][f][o.name] =
            @benchmarkable collision($p, $o)
    end
end

for (f, p) in zip(["straight", "magnetic"], particles)
    for (bname, bil) in zip(["mushroom", "psinai"], (bt, bt2))
        SUITE["next_collision"][f][bname] =
            @benchmarkable next_collision($p, $bil)
    end
end

for (f, p) in zip(["straight", "magnetic"], particles)
    for (bname, bil) in zip(["mushroom", "psinai"], (bt, bt2))
        ploc = copy(p) #location mutation screws up tests for different bts
        SUITE["bounce!"][f][bname] =
            @benchmarkable bounce!($ploc, $bil)
    end
end

let (f, p) = ("straight", particles[1])
    #resolvecollision! is indepent of particle type
    for o in obstacles
        ploc = copy(p)
        SUITE["resolvecollision!"][f][o.name] =
            @benchmarkable resolvecollision!($ploc, $o)
    end
end

for (f, p) in zip(["straight", "magnetic"], particles)
    ploc = copy(p)
    SUITE["propagate!"][f] =
        @benchmarkable propagate!($ploc, $proptime)
end

for (f, p) in zip(["straight", "magnetic"], particles)
    for (bname, bil) in zip(["mushroom", "psinai"], (bt, bt2))
        SUITE["lyapunovspectrum"][f][bname] =
            @benchmarkable lyapunovspectrum($p, $bil, 10000.0)
    end
end

################################################################################
## randominside – rewrite to add more non-collision benchmarks
################################################################################

SUITE["randominside"] = BenchmarkGroup() #I don't know any tags for randominside...
SUITE["randominside"]["straight"] = BenchmarkGroup(["straight"])

let f = "straight"
    for (bname, bil) in zip(["mushroom", "psinai"], (bt, bt2))
        SUITE["randominside"][f][bname] = @benchmarkable randominside($bil)
    end
end
