using Test
using DynamicalBilliards
using DynamicalBilliards.Testing

Floats = [Float16, Float32, Float64, BigFloat]

function alltypes(f, args...)
    r = 0.25; x = y = 1.0
    for T ∈ Floats
        bd = billiard_sinai(T(r), T(x), T(y))
        @test eltype(bd) == T

        p = randominside(bd)
        @test eltype(p) == T
        f(p, bd, args...)

        p = randominside(bd, T(1.0))
        @test eltype(p) == T
        f(p, bd, args...)
    end
end

function type_stability(p, bd)
    @testset "$(tag(p))" begin
        T = eltype(bd)
        for o in bd
            t, cp = collision(p, o)
            @test typeof(t) == T
            @test eltype(cp) == T
        end
        i, t, cp = next_collision(p, bd)
        @test typeof(t) == T
        @test eltype(cp) == T

        ct, poss, vels = evolve(p, bd, 5)
        @test eltype(ct) == T
        @test eltype(poss[1]) == eltype(vels[1]) == T
        ct, poss, vels = evolve(p, bd, 5.0)
        @test eltype(ct) == T
        @test eltype(poss[1]) == eltype(vels[1]) == T

        xt, yt, vxt, vyt, tt = timeseries(p, bd, T(10.0); dt = T(0.1))

        for x in (xt, yt, vxt, vyt, tt)
            @test eltype(x) == T
        end

    end
end

billiards_testset("type stability", type_stability; caller = alltypes)
