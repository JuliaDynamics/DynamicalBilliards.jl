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


# billiards_testset("Finite Collision time", test_finite_collisions)
# billiards_testset("min/max distance", test_maximum_distance, 1e-12)


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
