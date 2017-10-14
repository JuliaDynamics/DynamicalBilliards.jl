using Base.Test
using DynamicalBilliards

function type_stability(partnum=500; printinfo=true)
tim = time()
@testset "Type Stability: Random Sinai" begin
    Floats = [Float16, Float32, Float64, BigFloat]
    r = 0.25; x = y = 1.0
    @testset "Type: $(T), ω = $(ω)" for T ∈ Floats, ω ∈ [0, 1]
        bt = billiard_sinai(T(r), T(x), T(y), setting = "random")
        for obst in bt
            @test eltype(obst) == T
        end
        p = randominside(bt, T(ω))
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
        p = randominside(bt, T(ω))
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
    @testset "Straight Escape Time" begin
        bt = DynamicalBilliards.billiard_square_mushroom()
        for i in 1:partnum
            p = randominside(bt)

            et = escapetime(p, bt)
            @test et < Inf
        end#particle loop
    end
    @testset "Straight Escape Time" begin
        bt = DynamicalBilliards.billiard_square_mushroom()
        for i in 1:partnum
            p = randominside(bt, 2.0)

            et = escapetime(p, bt)
            @test et < Inf
        end#particle loop
    end
    if printinfo
        println("Results:")
        println("+ escapetime works for Particle and MagneticParticle")
        println("  and understands Doors.")
        println("+ The escape time is always finite.")
        println("+ randominside() work for FiniteWall and gives i.c.")
        println("+ inside mushroom cap.")
        println("+ Required time: $(round(time()-tim, 3)) sec.")
    end
    return
end
