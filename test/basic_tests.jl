using Test
using DynamicalBilliards
using DynamicalBilliards.Testing

function test_finite_collisions(p, bd, N = 1e6)
    n = 0
    i, t = next_collision(p, bd)
    s = Set(i)
    @testset "$(tag(p)) 0 < t < Inf" begin
        while n < N
            i, t = bounce!(p, bd)
            if 0 < t < Inf
                @test true
            else
                error("$(tag(p, bd)) collision time $(t) with obst.index $(i)")
            end
            n += 1
            push!(s, i)
        end
    end
    if typeof(p) <: Particle && !isperiodic(bd)
        iprev = 0
        n = 0
        @testset "collide with each obstacle ONCE!" begin
            while n < 1000
                i, t = bounce!(p, bd)
                if i != 4
                    if i != iprev
                        @test true
                    else
                        error("$(tag(p)) consecutive collisions with $(i)")
                    end
                end
                iprev = i
                n += 1
            end
        end
    end
    @testset "$(tag(p, bd)) collided with all obstacles" begin
        @test s == Set(1:length(bd))
    end
end

function test_maximum_distance(p, bd, thres = 1e-12, N = 1e6)
    acc = acclevel(p, bd)
    n = 0
    t, i = next_collision(p, bd)
    maxds = fill(-Inf, length(bd))
    minds = fill(Inf, length(bd))
    thresholds = fill(0, length(bd))
    isp = isperiodic(bd)
    @testset "$(tag(p, bd)) max/min distance" begin
        while n < N

            i, t, cp = next_collision(p, bd)
            d = distance(cp, bd[i])

            d < minds[i] && (minds[i] = d)
            d > maxds[i] && (maxds[i] = d)

            i, t = bounce!(p, bd)

            @test abs(d) < acc
            abs(d) > thres && (thresholds[i] += 1)

            n += 1
        end
        println(tag(p, bd))
        println("minds: ", minds)
        println("maxds: ", maxds)
        println("more than $(thres): ", thresholds, " (acc = $(acc))")
        println()
    end
end

billiards_testset("Finite Collision time", test_finite_collisions)
billiards_testset("min/max distance", test_maximum_distance, 1e-12)
