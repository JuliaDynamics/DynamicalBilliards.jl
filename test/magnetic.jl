using Base.Test
using DynamicalBilliards

function magnetic_sinai(partnum=500; printinfo = true)
tim = time()
@testset "Magnetic Sinai" begin
    for ω in [-0.5, 1.2]
        @testset "params: ω=$(ω), x=$(x), y=$(y), r=$(r)" for (r, x, y) in [(0.4, 1.0, 1.0), (0.3, 1.5, 1.0)]
            bt = billiard_sinai(r, x, y)
            c = bt[5].c
            tt=1000.0
            for i in 1:partnum
                p = randominside(ω, bt)

                t, poss, vels,  = evolve!(p, bt, tt)
                @test t[end] != Inf

                xt = [pos[1] for pos in poss]; yt = [pos[2] for pos in poss]
                mind = minimum(
                sqrt(((xt[i] - c[1])^2 + (yt[i] - c[2])^2)) for i in 1:length(xt))

                @test mind - r > -1e-15
                @test maximum(xt) ≤ x
                @test minimum(xt) ≥ 0
                @test maximum(yt) ≤ y
                @test minimum(yt) ≥ 0
            end#particle loop
        end#radius loop
    end#omega loop
end#testset
if printinfo
    println("Results:")
    println("+ evolve!() works for MagneticParticle.")
    println("+ randominside(ω) returns MagneticParticle.")
    println("+ relocate(), collisiontime(), resolvecollision() work for")
    println("  MagneticParticle with InfiniteWall and Disk.")
    println("+ particle never leaks the billiard table.")
    println("+ Required time: $(round(time()-tim, 3)) sec.")
end
return
end

function magnetic_periodic(partnum = 500; printinfo = true)
tim = time()
@testset "Magnetic Periodic Sinai" begin
    for ω in [0.2, -0.15]
        @testset "params: ω=$(ω), x=$(x), y=$(y), r=$(r)" for (r, x, y) in [(0.4, 1.0, 1.0), (0.42, 1.0, 1.2)]
            bt = billiard_sinai(r, x, y; setting="periodic")
            xmin, ymin, xmax, ymax = cellsize(bt)
            d = bt[5]
            c = d.c
            tt=1000.0
            invalid = 0
            minddist = min(x, y)

            for i in 1:partnum
                p = randominside(ω, bt)
                ts, poss, vels = evolve!(p, bt, tt)

                @test ts[end] != Inf

                mincolt = minimum(ts[3:end])
                error_level = 1e-8 #this huge error comes from the modulo operation
                xt = [mod(pos[1], xmax) for pos in poss]
                yt = [mod(pos[2], ymax) for pos in poss]

                mind = minimum(
                sqrt(((xt[i] - c[1])^2 + (yt[i] - c[2])^2)) for i in 1:length(xt))

                @test mind - d.r ≥ -error_level
                @test mincolt ≥ minddist-2r
            end#particle loop
        end#x,y loop (testset)
    end#omega loop
end#testset
@testset "Magnetic Periodic BigFloat" begin
    ω = big(0.2)
    (r, x, y) = big.([0.4, 1.0, 1.0])
    bt = billiard_sinai(r, x, y; setting="periodic")
    xmin, ymin, xmax, ymax = cellsize(bt)
    d = bt[5]
    c = d.c
    tt=1000.0
    invalid = 0
    minddist = min(x, y)

    for i in 1:1
        p = randominside(ω, bt)
        ts, poss, vels = evolve!(p, bt, tt)
        @test eltype(poss[1]) == BigFloat

        @test ts[end] != Inf

        mincolt = minimum(ts[3:end])
        error_level = 1e-8 #this huge error comes from the modulo operation
        xt = [mod(pos[1], xmax) for pos in poss]
        yt = [mod(pos[2], ymax) for pos in poss]

        mind = minimum(
        sqrt(((xt[i] - c[1])^2 + (yt[i] - c[2])^2)) for i in 1:length(xt))

        @test mind - d.r ≥ -error_level
        @test mincolt ≥ minddist-2r
    end#particle loop
end#testset
if printinfo
    println("Results:")
    println("+ billiard_sinai(setting=\"periodic\") and randominside() work.")
    println("+ relocate(), collisiontime(), resolvecollision() work for")
    println("  MagneticParticle with PeriodicWall.")
    println("+ collisiontime() ≤ min(x,y) - 2r.")
    println("+ Particle never invades the Disk.")
    println("+ Collision time is never Infinite (no pinned particles).")
    println("+ All the above also work for BigFloat.")
    println("+ Required time: $(round(time()-tim, 3)) sec.")
end
return
end#function
