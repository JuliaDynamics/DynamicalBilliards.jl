using Test
using DynamicalBilliards
using DynamicalBilliards.Testing

function test_no_escape(p, bd, N = 1e4)
    xmin, ymin, xmax, ymax = cellsize(bd)

    ct, poss, vels = evolve!(p, bd, Int(N))
    @test length(ct) == N + 1
    @test typeof(poss) == typeof(vels) == Vector{SVector{2, eltype(p)}}

    if typeof(p) <: MagneticParticle
        xt, yt = construct(ct, poss, vels, p.ω, 0.1)
    else
        xt, yt = construct(ct, poss, vels)
    end

    @test typeof(xt) == typeof(yt) == Vector{eltype(p)}

    dx1 = minimum(xt) - xmin
    dx2 = xmax - maximum(xt)
    dy1 = minimum(yt) - ymin
    dy2 = ymax - maximum(yt)

    println(tag(p, bd))
    println("dx1: ", dx1, " dx2: ", dx2)
    println("dy1: ", dy1, " dy2: ", dy2)
    println()
    @testset "$(tag(p, bd)) no escape" begin
        for d in (dx1, dx2, dy1, dy2)
            @test d ≥ -acclevel(p, bd)
        end
    end
end

billiards_testset("No escape", test_no_escape; caller = ergodic_tests)



function test_movin_periodic(p, bd, N = 1e2)
    xmin, ymin, xmax, ymax = cellsize(bd)

    ct, poss, vels = evolve!(p, bd, N)
    @test length(ct) != N + 1
    @test typeof(poss) == typeof(vels) == Vector{SVector{2, eltype(p)}}

    if typeof(p) <: MagneticParticle
        xt, yt = construct(ct, poss, vels, p.ω, 0.1)
    else
        xt, yt = construct(ct, poss, vels)
    end

    @test typeof(xt) == typeof(yt) == Vector{eltype(p)}

    dx1 = minimum(xt) - xmin
    dx2 = xmax - maximum(xt)
    dy1 = minimum(yt) - ymin
    dy2 = ymax - maximum(yt)

    @test maximum(abs.((dx1, dx2))) > xmax - xmin
    @test maximum(abs.((dy1, dy2))) > ymax - ymin

end

billiards_testset("movin periodic", test_movin_periodic; caller = periodic_tests)
