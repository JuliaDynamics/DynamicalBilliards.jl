using Test
using DynamicalBilliards
using DynamicalBilliards.Testing

function coordinate_invariance(p, bd, args...)
    @testset "$(tag(p, bd)) Coordinate change" begin
        for t in 1:2e4
            i, tmin = bounce!(p, bd)
            ξ, sφ = to_bcoords(p, bd[i])
            pos, vel = from_bcoords(ξ, sφ, bd[i])
            @test *(isapprox.(p.pos, pos, atol=1e-4)...)
            @test *(isapprox.(p.vel, vel, atol=1e-4)...)
        end
    end
end

billiards_testset("Coordinate change invariance", coordinate_invariance; caller = omni_tests_noellipse)



function cut_psos(args...)
    t = 1000
    bt = billiard_sinai(0.25)
    plane = InfiniteWall([0.5, 0.0], [0.5, 1.0], [-1.0, 0.0])
    @testset "PSOS: pinned particle" begin

        # case of pinned that crosses plane but not walls
        p = MagneticParticle(0.2, 0.5, -π/2, 1/0.3)
        a, b = psos(bt, plane, t, p)
        @test typeof(a) <: Vector{<:SVector}
        @test length(a) == length(b) == 1

        # case of pinned that doesn't cross anything
        p = MagneticParticle(0.1, 0.5, -π/2, 1/0.05)
        a, b = psos(bt,plane, t, p)
        @test typeof(a) <: Vector{<:SVector}
        @test length(a) == length(b) == 0

        # case of pinned but does cross periodic walls *and* section
        bd = billiard_sinai(;setting = "periodic")
        p = MagneticParticle(0.25, 0.25, -π/2 + π/4, 0.44*2)
        a, b = psos(bd, plane, t, p)
        @test typeof(a) <: Vector{<:SVector}
        @test length(a) == length(b) == 1
    end

    for ω ∈ [0.0, 1.0]
        p = ω == 0 ? randominside(bt) : randominside(bt, ω)
        @testset "$(tag(p, bt)) psos" begin
            a, b = psos(bt, plane, t, p)
            for j in 1:length(a)
                @test a[j][1] ≈ 0.5
                @test 0 < a[j][2] < 1
                @test -1 < b[j][1] < 1
                @test -1 < b[j][2] < 1
            end
        end
    end

    @testset "many particles" begin
        a, b = psos(bt, plane, 1000, 10)
        for i in 1:length(a)
            @test length(a[i]) == length(b[i])
        end
    end
end

billiards_testset("PSOS", identity; caller = cut_psos)




function fills_boundarymap(p, bd)
    @testset "$(tag(p, bd)) Fills boundary map" begin
        bmap, = boundarymap(p, bd, 10000)
        ξs = [b[1] for b in bmap]; sφs = [b[2] for b in bmap]
        partition = (10,10) # partition size
        A = falses(partition)
        l = totallength(bd)
        ε = (l/partition[1], 2/partition[2]) # box size
        c = 0
        for point ∈ zip(ξs, sφs)
            id = clamp.(ceil.(Int, (point .- (0, -1))./ε), (1,1), partition)
            if !A[id...]
                A[id...] = true
                c += 1
                if c == partition[1]*partition[2]
                    break
                end
            end
        end
        for a in A
            @test a
        end
    end
end

billiards_testset("Fills boundary map", fills_boundarymap; caller = ergodic_tests)


function phasespace_ratio(f, g)
    @testset "Bunimovich" begin
        t = 100000.0
        bt = billiard_bunimovich()
        for ω in [0.0, 0.1]
            p = ω == 0 ? randominside(bt) : randominside(bt, ω)
            @testset "$(tag(p, bt))" begin
                φ = π/4 * rand() # so that we never find bouncing walls
                p.vel = (cos(φ), sin(φ))
                ratio, = f(bt,t, randominside(bt), 0.1)
                @test ratio == 1.0
            end
        end
    end

    t = 1000000.0
    l = 1.0; r = 1.0
    @testset "Mushroom w=$(w)" for w ∈ [0.2, 0.4]

        bt = billiard_mushroom(l, w, r)
        @testset "regular" begin
            p = MushroomTools.randomregular(l, w, r)
            ratio, = f(bt,t, p, 0.1)
            trueratio =  1 - g(l,w,r)
            # Only one regular particle covers very small amount of space:
            @test ratio < trueratio
        end
        @testset "chaotic" begin
            p = MushroomTools.randomchaotic(l, w, r)
            ratio, = f(bt, t, p, 0.1)
            trueratio = g(l,w,r)
            @test trueratio - 0.1 ≤ ratio ≤ trueratio + 0.1
        end
    end
end
ratio_caller(x, f, g) = phasespace_ratio(f, g)

billiards_testset("2D Phasespace Ratio", identity, boundarymap_portion, MushroomTools.g_c_2D;
caller = ratio_caller)

billiards_testset("3D Phasespace Ratio", identity, phasespace_portion, MushroomTools.g_c_3D;
caller = ratio_caller)
