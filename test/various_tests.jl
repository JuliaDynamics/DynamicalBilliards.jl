using Test
using DynamicalBilliards
using DynamicalBilliards.Testing

function escapetest(args...)
    l = 1.0; w = 0.2; r = 1.0
    bd = billiard_mushroom(l, w, r)

    for i in 1:20
        p = MushroomTools.randomchaotic(l, w, r)

        t = escapetime(p, bd, Int(1e6))
        @test t < Inf

        p = MushroomTools.randomregular(l, w, r)
        t = escapetime(p, bd, Int(1e6))
        @test t == Inf
    end
end

billiards_testset("escape time", identity; caller = escapetest)

function meancoltest(args...)
    p = Particle(0.1, 0.1, 0.1)
    @testset "Sinai" begin
        for r in [0.25, 0.35]
            bd = billiard_sinai(r)
            a = π*(1 - π*(r)^2)/(4 + 2π*r)
            m = meancollisiontime(p, bd, Int(1e6))
            @test a ≈ m rtol = 0.1
        end
    end
    @testset "Stadium" begin
        l = 1.0; w = 1.0
        bd = billiard_stadium(l, w)
        a = π*(w*l + π*(w/2)^2)/(2l + π*w)
        m = meancollisiontime(p, bd, Int(1e6))
        @test a ≈ m rtol = 0.1
    end
    @testset "Stadium" begin
        l = 1.0; w = 1.0
        bd = billiard_stadium(l, w)
        p = MagneticParticle(0.1, 0.1, 0.1, 1.0)
        a = π*(w*l + π*(w/2)^2)/(2l + π*w)
        m = meancollisiontime(p, bd, Int(1e6))
        @test a ≈ m rtol = 0.1
    end
end

billiards_testset("mean collision time", identity; caller = meancoltest)

function noparticle_interafaces(args...)
    bd = billiard_mushroom()
    @test typeof(escapetime(bd, 100)) == Float64
    @test typeof(meancollisiontime(bd, 100)) == Float64
    @test typeof(evolve(bd, 100)[1][1]) == Float64
    # @test typeof(lyapunovspectrum(bd, 100.0)[1]) == Float64
end

billiards_testset("no-particle interface", identity; caller = noparticle_interafaces)


function ispinned_tests(args...)
    bd = billiard_sinai(;setting = "periodic")
    # case of pinned that crosses plane but not walls
    p = MagneticParticle(0.2, 0.5, -π/2, 1/0.3)
    @test ispinned(p, bd)
    # case of pinned that doesn't cross anything
    p = MagneticParticle(0.1, 0.5, -π/2, 1/0.05)
    @test ispinned(p, bd)
    # case of pinned but does cross periodic walls *and* section
    p = MagneticParticle(0.0, 1.0, -π/2, 0.44*2)
    @test ispinned(p, bd)
    # Sanity check for 2 nonpinned particles
    @test !ispinned(randominside(bd), bd)
    p = MagneticParticle(0.2, 0.8, -π/2, 1/0.3)
    @test !ispinned(p, bd)
end

billiards_testset("ispinned", identity; caller = ispinned_tests)

# Quick an sloppy parallelization test
function test_parallelized(args...)
    bd = billiard_mushroom()
    N = 100
    particles = [randominside(bd) for i in 1:N]
    for f in (meancollisiontime!, escapetime!, lyapunovspectrum!)
        res = parallelize(f, bd, N, particles)
        @test length(res) == N
    end

    res, ai = parallelize(boundarymap, bd, N, particles)
    @test length(res) == N
end

billiards_testset("parallelize", identity; caller = test_parallelized)
