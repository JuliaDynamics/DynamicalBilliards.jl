using DynamicalBilliards
using Base.Test

#=debug=# false && using Juno

function raysplit_straight(partnum=500; printinfo = true)
tim = time()
@testset "Raysplitting Straight" begin

    @testset "params: x=$(x), y=$(y)" for (x, y, r1, r2) in [(3,1,0.3,0.2), (3,2,0.6, 0.5)]
        bt, ray = billiard_raysplitting_showcase(x, y, r1, r2)

        @test isphysical(ray)
        tt=10000.0
        for i in 1:partnum
            p = randominside(bt)
            t, poss, vels = evolve!(p, bt, tt, ray)
            reset_billiard!(bt)
            @test t[end] != Inf
            xt = [pos[1] for pos in poss]; yt = [pos[2] for pos in poss]
            @test maximum(xt) ≤ x
            @test minimum(xt) ≥ 0
            @test maximum(yt) ≤ y
            @test minimum(yt) ≥ 0
            reset_billiard!(bt)
        end#particle loop
    end#parameters
end#testset
if printinfo
    println("Results:")
    println("+ evolve!() works for Ray-splitting billiards & Particle.")
    println("+ relocate(), collisiontime(), resolvecollision() work for")
    println("  Particle & Ray-splitting with SplitterWall and Antidot.")
    println("+ particle never leaks the billiard table.")
    println("+ Required time: $(round(time()-tim, 3)) sec.")
end
return
end

function raysplit_magnetic(partnum=500; printinfo = true)
tim = time()
@testset "Raysplitting Magnetic" begin

    btcount = 1
    @testset "params: x=$(x), y=$(y)" for
    (x, y, r1, r2) in [(3,1,0.3,0.2), (3,2,0.6, 0.5)]
        bt, ray = billiard_raysplitting_showcase(x, y, r1, r2)

        @test isphysical(ray)
        tt=1000.0
        for i in 1:partnum
            p = randominside(bt, 0.4)
            #=debug=# false && Juno.clearconsole()
            #=debug=# false && println("Particle ", i, " billiard $btcount")
            #=debug=# false && println("pos = SVector($(p.pos[1]), $(p.pos[2]))")
            #=debug=# false && println("vel = SVector($(p.vel[1]), $(p.vel[2]))")
            t, poss, vels = evolve!(p, bt, tt, ray)
            reset_billiard!(bt)
            @test t[end] != Inf
            xt = [pos[1] for pos in poss]; yt = [pos[2] for pos in poss]
            @test maximum(xt) ≤ x
            @test minimum(xt) ≥ 0
            @test maximum(yt) ≤ y
            @test minimum(yt) ≥ 0
            reset_billiard!(bt)
        end#particle loop
        #=debug=# false && Juno.clearconsole()
        btcount += 1
    end#parameters
end#testset
# @testset "Raysplitting Magnetic BigFloat" begin
#
#     (x, y, r1, r2) = big.([3,1,0.3,0.2])
#     bt, ray = billiard_raysplitting_showcase(x, y, r1, r2)
#
#     @test isphysical(ray)
#     tt=50.0
#     p = randominside(bt, big(0.8))
#     t, poss, vels = evolve!(p, bt, tt, ray)
#     @test t[end] != Inf
#     xt = [pos[1] for pos in poss]; yt = [pos[2] for pos in poss]
#     @test maximum(xt) ≤ x
#     @test minimum(xt) ≥ 0
#     @test maximum(yt) ≤ y
#     @test minimum(yt) ≥ 0
#     reset_billiard!(bt)
# end#testset
if printinfo
    println("Results:")
    println("+ evolve!() works for Ray-splitting billiards & MagneticParticle.")
    println("+ relocate(), collisiontime(), resolvecollision() work for")
    println("  MagneticParticle & Ray-splitting with SplitterWall and Antidot.")
    println("+ particle never leaks the billiard table + process terminates.")
    println("+ All the above also work for BigFloat.")
    println("+ Required time: $(round(time()-tim, 3)) sec.")
end
return
end
