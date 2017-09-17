using Base.Test
using DynamicalBilliards

function straight_sinai(partnum=500; printinfo=true)
tim = time()
@testset "Straight Sinai" begin
    @testset "rectangle params: $((r, x, y))" for (r, x, y) ∈ [(0.15, 1.0, 1.0), (0.3, 1.5, 1.0), (0.4, 1.4, 2.2)]
        bt = billiard_sinai(r, x, y)
        d = bt[5]
        c = d.c
        r = d.r
        tt = 10000.0

        for i in 1:partnum
            p = randominside(bt)
            ts, poss, vels = evolve!(p, bt, tt)

            error_level = 1e-15

            xt = [pos[1] for pos in poss]; yt = [pos[2] for pos in poss]
            dist = @. sqrt(((xt - c[1])^2 + (yt - c[2])^2))
            mind = minimum(dist)
            @test mind - r > -error_level

            @test maximum(xt) <= x
            @test minimum(xt) >= 0
            @test maximum(yt) <= y
            @test minimum(yt) >= 0
        end#particle loop
    end#x,y loop
end#testset
if printinfo
    println("Results:")
    println("+ evolve!() works for Particle and has finite collision time.")
    println("+ billiard_sinai() and randominside() work.")
    println("+ relocate(), collisiontime(), resolvecollision() work for")
    println("  for Particle with FiniteWall and Disk.")
    println("+ particle never leaks the billiard table.")
    println("+ Required time: $(round(time()-tim, 3)) sec.")
end
return
end#function


function straight_periodic(partnum = 500; printinfo = true)
tim = time()
@testset "Straight Periodic Sinai" begin
    @testset "rectangle dims $((r, x, y))" for (r, x, y) in [(0.3, 1.5, 1.0), (0.5, 1.4, 2.2), (0.2, 0.8, 2.2)]
        bt = billiard_sinai(r, x, y; setting="periodic")
        xmin, ymin, xmax, ymax = cellsize(bt)
        d = bt[5]
        c = d.c
        tt=10000.0
        invalid = 0
        minddist = min(x, y)

        for i in 1:partnum
            mincolt = Inf
            p = randominside(bt)
            ts, poss, vels = evolve!(p, bt, tt)

            if length(ts) > 2
                mincolt = minimum(ts[3:end])
            else
                continue
            end

            error_level = 1e-8 #this huge error comes from the modulo operation
            xt = [mod(pos[1], xmax) for pos in poss]
            yt = [mod(pos[2], ymax) for pos in poss]

            mind = minimum(
            sqrt(((xt[i] - c[1])^2 + (yt[i] - c[2])^2)) for i in 1:length(xt))

            @test mind - d.r ≥ -error_level
            @test mincolt ≥ minddist-2r

        end#particle loop
    end#x,y loop
end#testset
if printinfo
    println("Results:")
    println("+ evolve!() works for Particle and has finite collision time.")
    println("+ billiard_sinai(setting=\"periodic\") and randominside() work.")
    println("+ relocate(), collisiontime(), resolvecollision() work for")
    println("  for Particle with PeriodicWall and Disk.")
    println("+ collisiontime() ≤ min(x,y) - 2r.")
    println("+ Particle never invades the Disk.")
    println("+ Required time: $(round(time()-tim, 3)) sec.")
end
return
end#function
