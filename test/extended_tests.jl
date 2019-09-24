using Test
using DynamicalBilliards
using DynamicalBilliards.Testing

function test_no_escape(p, bd, N = 1e4)
    xmin, ymin, xmax, ymax = cellsize(bd)

    xt, yt, vyt, vxt, tvec = timeseries(p, bd, Int(N); dt = Inf)

    @test length(xt) == N + 1
    @test typeof(xt) == typeof(tvec) == Vector{eltype(p)}

    @testset "$(tag(p, bd)) correct tvector" begin
        for i in 2:length(tvec)
            @test tvec[i] > tvec[i-1]
        end
    end

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

function test_visited_obstacles(p, bd, N = 50)

    ts, obst = visited_obstacles(p, bd, N)
    cts, = evolve(p, bd, N)

    @test length(ts) == length(cts)
    cc = cumsum(cts)
    for i in eachindex(cc)
        @test cc[i] ≈ ts[i]  atol = 1e-12
    end
end
billiards_testset("Visited obstacles", test_visited_obstacles; caller = ergodic_tests)



function test_movin_periodic(p, bd, N = 1e3)
    xmin, ymin, xmax, ymax = cellsize(bd)

    xt, yt = timeseries(p, bd, Int(N), dt = Inf)
    @test length(xt) == N + 1
    @test typeof(xt) == typeof(xt) == Vector{eltype(p)}

    dx1 = minimum(xt) - xmin
    dx2 = xmax - maximum(xt)
    dy1 = minimum(yt) - ymin
    dy2 = ymax - maximum(yt)

    @test maximum(abs.((dx1, dx2))) > xmax - xmin
    @test maximum(abs.((dy1, dy2))) > ymax - ymin

end

billiards_testset("movin periodic", test_movin_periodic; caller = periodic_tests)

function test_raysplit_ts_inside(p, bd, ray, N = 1e4)

    xmin, ymin, xmax, ymax = cellsize(bd)

    xt, yt = timeseries(p, bd, Int(N), ray, dt = Inf)
    @test length(xt) == N + 1
    @test typeof(xt) == typeof(xt) == Vector{eltype(p)}

    dx1 = minimum(xt) - xmin
    dx2 = xmax - maximum(xt)
    dy1 = minimum(yt) - ymin
    dy2 = ymax - maximum(yt)

    for d in (dx1, dx2, dy1, dy2)
        @test d ≥ -acclevel(p, bd)
    end

end

billiards_testset("timeseries raysplit", test_raysplit_ts_inside; caller = ray_tests)
