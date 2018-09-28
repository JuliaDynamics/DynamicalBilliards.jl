using Test
using DynamicalBilliards
using DynamicalBilliards.Testing

"""
    test_finite_collisions(p, bd, N = 1e4)

1. All `bounce!` calls result in finite, non-zero collision time.
1. `bounce!` is tested with all non-periodic obstacle types.
2. Particles collides with all obstacles in `bd`.
3. For `Particle` tests that it collides once with obstacles (other than `Semicircle`).
"""
function test_finite_collisions(p, bd, N = 1e4)
    n = 0
    t, i = next_collision(p, bd)
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

"""
    test_maximum_distance(p, bd, level = 1e-9, N = 1e4)

1. Measure max/min distances of particle with collided obstacle.
2. Test that each of them does not exceed `level` (in magnitude).
3. Print the distances for reference.
"""
# function test_maximum_distance(p, bd, level = 1e-3, thres = 1e-9, N = 1e6)
#     n = 0
#     t, i = next_collision(p, bd)
#     maxds = fill(-Inf, length(bd))
#     minds = fill(Inf, length(bd))
#     thresholds = fill(0, length(bd))
#     isp = isperiodic(bd)
#     @testset "$(tag(p, bd)) max/min distance" begin
#         while n < N
#             i, t = bounce!(p, bd)
#             d = distance(p, bd[i])
#             d < minds[i] && (minds[i] = d)
#             d > maxds[i] && (maxds[i] = d)
#             if !isperiodic(bd) || i ∉ bd.peridx
#                 @test abs(d) < level
#                 abs(d) > thres && (thresholds[i] += 1)
#             end
#             n += 1
#         end
#     end
#     println(tag(p, bd))
#     println("minds: ", minds)
#     println("maxds: ", maxds)
#     println("distances more than $(thres): ", thresholds)
#     println()
# end

function test_maximum_distance(p, bd, level = 1e-8, thres = 1e-9, N = 1e6)
    n = 0
    t, i = next_collision(p, bd)
    maxds = fill(-Inf, length(bd))
    minds = fill(Inf, length(bd))
    thresholds = fill(0, length(bd))
    isp = isperiodic(bd)
    @testset "$(tag(p, bd)) max/min distance" begin
        while n < N

            t, i = next_collision(p, bd)
            pos = propagate_pos(p.pos, p, t)
            d = distance(pos, bd[i])

            d < minds[i] && (minds[i] = d)
            d > maxds[i] && (maxds[i] = d)

            i, t = bounce!(p, bd)

            if !isperiodic(bd) || i ∉ bd.peridx
                @test abs(d) < level
                abs(d) > thres && (thresholds[i] += 1)
            end
            n += 1
        end
    end
    println(tag(p, bd))
    println("minds: ", minds)
    println("maxds: ", maxds)
    println("distances more than $(thres): ", thresholds)
    println()
end



billiards_testset("Finite Collision time", test_finite_collisions, 1e6)
billiards_testset("min/max distance", test_maximum_distance, 1e-8, 1e-9)


# STEP 1 : Add maximum distance AFTER JUST CALLING next_collision and
# propagate! What is the maximum distance just from these calls? between
# particle and point of collision?
# Should I transform all functions to return time to collisoin as well as
# Collision point? Then I can move particle to collision point and then
# relocate with respect to the normal.

#
# 
# function test_maximum_relocations(p, bd, upper = 5, N = 1e6)
#     n = 0
#     maxk = zeros(Int, length(bd))
#     while n < N
#         tmin::Float64, i::Int = next_collision(p, bd)
#         if tmin == -Inf
#             error("tmin = $(tmin) for i = $(i)")
#         end
#         o = bd[i]
#         tmin, k = relocate!(p, o, tmin)
#         resolvecollision!(p, o)
#         c = Int(log2(k))
#         maxk[i] < c && (maxk[i] = c)
#         n += 1
#     end
#     println("$(tag(p,bd)) max relocation loops:")
#     println(maxk)
#     println()
# end
#
# p, mp = testparticles()
# bd = omnibilliard()
# test_maximum_relocations(mp, bd)

# billiards_testset("max relocation loops", test_maximum_relocations, 5)

# semicircle seems terrible for the relocation algorithm...?
# does 45 recursive steps wtf?
# I have to rework what the fuck is going on inside there.

# Put in a limiter k < 5 and after that do a special relocation
# that moves along the normal vector at collision ?
# Yeap!


# This seems like a really good plan!
