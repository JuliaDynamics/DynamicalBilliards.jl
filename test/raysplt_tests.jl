using Test
using DynamicalBilliards
using DynamicalBilliards.Testing

# %% Basic tests

function simple_raysplit(f, args...)
    bd, ray = basic_ray(true)
    p = randominside(bd)
    f(p, bd, ray, args...)

    bd, ray = basic_ray(false)
    mp = randominside(bd, 2.0)
    f(mp, bd, ray, args...)
end


function test_finite_collisions_ray(p, bd, ray, N = 5e5)
    n = 0
    i, t = next_collision(p, bd)
    s = Set(i)
    raysidx = DynamicalBilliards.raysplit_indices(bd, ray)

    @testset "$(tag(p)) 0 < t < Inf" begin
        while n < N
            i, t = bounce!(p, bd, raysidx, ray)
            if 0 < t < Inf
                @test true
            else
                error("$(tag(p, bd, ray)) collision time $(t) with obst.index $(i)")
            end
            n += 1
            push!(s, i)
        end
    end
    @testset "$(tag(p, bd)) collided with all obstacles" begin
        @test s == Set(1:length(bd))
    end
end


function ray_maximum_distance(p, bd, ray, thres = 1e-13, N = 5e5)
    acc = acclevel(p, bd, ray)
    n = 0
    t, i = next_collision(p, bd)
    maxds = fill(-Inf, length(bd))
    minds = fill(Inf, length(bd))
    thresholds = fill(0, length(bd))
    isp = isperiodic(bd)
    raysidx = DynamicalBilliards.raysplit_indices(bd, ray)

    @testset "$(tag(p, bd, ray)) max/min distance" begin
        while n < N

            i, t, cp = next_collision(p, bd)
            d = distance(cp, bd[i])

            d < minds[i] && (minds[i] = d)
            d > maxds[i] && (maxds[i] = d)

            i, t = bounce!(p, bd, raysidx, ray)

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

billiards_testset("Ray Collision time", test_finite_collisions_ray; caller = simple_raysplit)
billiards_testset("Ray min/max distance", ray_maximum_distance; caller = simple_raysplit)



# %% Extreme angles test

function extreme_raysplit(f, args...)
    bd, ray = extreme_ray(true)
    p = randominside(bd)
    f(p, bd, ray, args...)

    bd, ray = extreme_ray(false)
    mp = randominside(bd, 2.0)
    f(mp, bd, ray, args...)
end

billiards_testset("ERay Collision time", test_finite_collisions_ray; caller = extreme_raysplit)
billiards_testset("ERay min/max distance", ray_maximum_distance; caller = extreme_raysplit)


# %% Totally brutal tests

function inside_antidot(args...)
    bd, ray = extreme_ray(true)
    pa = randominside(bd)
    mp = randominside(bd, 1.0)
    isray(i) = (i == 1 || i == 2)
    raysidx = DynamicalBilliards.raysplit_indices(bd, ray)

    for p in (pa, mp); @testset "$(tag(p, bd, ray)) BRUTAL" begin
        bd, ray = extreme_ray(typeof(p) <: Particle)
        n = 0
        iprev = 0
        check = false

        while n < 5e5
            i, t = bounce!(p, bd, raysidx, ray)

            if check
                @test i == iprev
            end

            check = isray(i) && !(bd[i].pflag)
            iprev = i

            n += 1
        end
    end; end
end

billiards_testset("RAY Stays inside", identity; caller = inside_antidot)
