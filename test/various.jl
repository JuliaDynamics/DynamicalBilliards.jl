using Test
using DynamicalBilliards

function type_stability(partnum=500; printinfo=true)
tim = time()
@testset "Type Stability: Random Sinai" begin
    Floats = [Float16, Float32, Float64, BigFloat]
    r = 0.25; x = y = 1.0
    @testset "Type: $(T), ω = $(ω)" for T ∈ Floats, ω ∈ [0.0, 1.0]
        bt = billiard_sinai(T(r), T(x), T(y), setting = "random")
        for obst in bt
            @test eltype(obst) == T
        end
        p = ω == 0 ? randominside(bt) : randominside(bt, T(ω))
        if ω == 0
            @test typeof(p) <: Particle
        else
            @test typeof(p) <: MagneticParticle
        end
        @test eltype(p) == T
        tt = T(100.0)
        ts, poss, vels = evolve!(p, bt, tt)
        @test eltype(ts) == T
        @test eltype(poss[1]) == T

    end#type loop
end#testset
@testset "Type Stability: Hexagon" begin
    Floats = [Float16, Float32, Float64, BigFloat]
    r = 0.75; x = y = 1.0
    @testset "Type: $(T), ω = $(ω)" for T ∈ Floats, ω ∈ [0, 0.2]
        bt = billiard_hexagonal_sinai(T(0.4), T(0.6), setting="periodic")
        for obst in bt
            @test eltype(obst) == T
        end
        p = ω == 0 ? randominside(bt) : randominside(bt, T(ω))
        if ω == 0
            @test typeof(p) <: Particle
        else
            @test typeof(p) <: MagneticParticle
        end
        @test eltype(p) == T
        tt = T(100.0)
        ts, poss, vels = evolve!(p, bt, tt)
        @test eltype(ts) == T
        @test eltype(poss[1]) == T

    end#type loop
end#testset
if printinfo
    println("Results:")
    println("+ billiard_sinai() and billiard_polygon() conserve types")
    println("+ randominside() returns correct types for ω=0 and ω=1")
    println("+ evolve!() conserves type for all the combinations of the above.")
    println("+ evolve!() works for Random sinai and Hexagonal periodic billiard.")
end
return
end#function

function escape_times(partnum=500; printinfo=true)
    tim = time()
    bt = billiard_mushroom()
    @testset "Straight Escape Time" begin
        port = 0
        for i in 1:partnum
            p = randominside(bt)

            et = escapetime(p, bt, 10000, warning=false)
            @test typeof(et)==Float64
            if et == Inf
                port +=1
            end
        end#particle loop
        @test port < partnum
    end
    @testset "Magnetic Escape Time" begin
        port = 0
        for i in 1:partnum
            p = randominside(bt, 0.1)

            et = escapetime(p, bt, 10000, warning=false)
            @test typeof(et)==Float64
            if et == Inf
                port +=1
            end
        end#particle loop
        @test port < partnum
    end
    if printinfo
        println("Results:")
        println("+ escapetime works for Particle and MagneticParticle")
        println("  and understands Doors.")
        println("+ The escape time is always finite.")
        println("+ randominside() works for billiard_mushroom()!")
        println("+ Required time: $(round(time()-tim, digits=3)) sec.")
    end
    return
end


function meancoltimes(partnum=500; printinfo=true)
    tim = time()
    @testset "Mean Collision Times" begin
    bt = billiard_sinai()
    @testset "Sinai ω = $ω" for ω in (0.0, 0.5)
        for i ∈ 1:partnum
            p = ω == 0 ? randominside(bt) : randominside(bt, ω)

            κ1 = meancollisiontime!(p, bt, 1000)
            @test κ1 < 1.5
        end
    end
    j = 1
    bt = billiard_sinai(0.15, setting = "periodic")
    @testset "Periodic sinai ω = $ω" for ω ∈ [0.01, 0.88, 2.0]
        mcts = [3.097297, 2.201, 1.76]
        port = 0
        for i in 1:partnum
            p = ω == 0 ? randominside(bt) : randominside(bt, ω)

            mct = meancollisiontime!(p, bt, 1000)
            if mct == Inf
                port +=1
            else
                @test mcts[j] - 2.0 ≤ mct ≤ mcts[j] + 4.0
            end
        end#particle loop
        @test port < partnum
        partnum > 10 && ω != 0.01 && @test port > 0
        j += 1
    end
    end
    if printinfo
        println("Results:")
        println("+ mean collision time works.")
        println("+ Required time: $(round(time()-tim, digits=3)) sec.")
    end
    return
end

function _coordinates(bd, partnum=500; tolerance = 1e-5)
    intervals = arcintervals(bd)
    @testset "real -> boundary -> real" begin
        for i ∈ 1:partnum
            #create real position & velocity on boundary
            p = randominside(bd);
            tmin, i = next_collision(p, bd)
            relocate!(p, bd[i], tmin)

            #transform to Boundary coordinates
            ξ = intervals[i] + to_bcoords(p, bd[i])
            sφ = sin(DynamicalBilliards.incidence_angle(p, bd[i]))
            
            #transform back to real space coordinates
            rpos, rvel, ri = from_bcoords(ξ, sφ, bd, return_obstacle=true,
                                              intervals = intervals)

            @test i == ri
            @test norm(rpos-p.pos) < tolerance
            @test norm(rvel-p.vel) < tolerance
        end
    end
    
    @testset "boundary -> real -> boundary" begin
        for i ∈ 1:partnum
            #create to_bcoords & angle of incidence
            ξ = rand()*intervals[end]
            sφ = 2*rand() - 1.0
            
            #transform to real coordinates
            pos, vel, i = from_bcoords(ξ, sφ, bd, return_obstacle=true,
                                           intervals = intervals)
            p = Particle(pos, vel, SVector(0.0,0.0))
            #transform back to Birkhoff coordinates
            rξ = to_bcoords(pos, bd[i]) + intervals[i]
            rsφ = sin(DynamicalBilliards.incidence_angle(p, bd[i]))

            @test rξ <= intervals[end] && rξ >= 0
            
            @test abs(rξ - ξ) < tolerance
            @test abs(rsφ - sφ) < tolerance
        end
    end
end

function coordinates(partnum=500; printinfo=true)
    tim = time()

    tolerance = 1e-10
    
    @testset "Mushroom" begin
        bd = billiard_mushroom()
        _coordinates(bd, partnum, tolerance = tolerance)
    end

    @testset "Stadium" begin
        bd = billiard_mushroom()
        _coordinates(bd, partnum, tolerance = tolerance)
    end

    if printinfo
        println(
            """Results:
+ from_bcoords and to_bcoords return sane 
  results on mushrooms and stadiums
+ from_bcoords and to_bcoords are inverse
  functions of each other
""")
    end
    return
end
    
