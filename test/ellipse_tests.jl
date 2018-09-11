using Test
using DynamicalBilliards

function ellipse_tests(partnum=500; printinfo = true)
tim = time()
@testset "Ellipse Arclength" begin
    o1 = Ellipse([0,0], 1.0, 0.5)
    o2 = Ellipse([0,0], 0.5, 1.0)

    for x in [o1, o2]
        for y in [o1, o2]
            @test ellipse_arclength(π/2, x) == y.arc
        end
        @test ellipse_arclength(π, x) == 2x.arc
        @test ellipse_arclength(2π, x) == 4x.arc
        @test ellipse_arclength(π/2, x) + 2x.arc == ellipse_arclength(3π/2, x)
    end

    phis = 0.0:0.01:2π

    for i in 2:length(phis)
        @test ellipse_arclength(phis[i], o1) > ellipse_arclength(phis[i-1], o1)
    end

    # Arclength should be same as circle!
    o3 = Ellipse([0,0], 1.0, 1.0)
    d = Disk([0,0], 1.0)

    # @test totallength(o3) ≈ totallength(d)

end#testset

@testset "Iris Standard" begin
    bd = billiard_iris()
    o = bd[1]
    for i in 1:partnum
        p = randominside(bd)
        xt, yt = construct(evolve(p, bd, 1000)...)

        @test minimum(distance(SVector(x, y), o) for (x,y) in zip(xt, yt)) ≥ 0
        @test maximum(xt) <= 1
        @test minimum(xt) >= 0
        @test maximum(yt) <= 1
        @test minimum(yt) >= 0
    end
end

@testset "Iris Std. Ray" begin
    bd = billiard_iris()
    o = bd[1]
    refraction(φ, pflag, ω) = pflag ? 0.5φ : 2.0φ
    transmission(φ, pflag, ω) = pflag ? 1.0 : abs(φ) < π/4 ? 1.0 : 0.0
    raya = RaySplitter([1], transmission, refraction)
    for i in 1:partnum
        p = randominside(bd)
        xt, yt = construct(evolve(p, bd, 1000, raya)...)
        @test maximum(xt) <= 1
        @test minimum(xt) >= 0
        @test maximum(yt) <= 1
        @test minimum(yt) >= 0
    end
end

@testset "Equiv to Sinai" begin
    bd2 = billiard_sinai(0.25)
    bd1 = billiard_iris(0.25, 0.25)

    for i in 1:10
        N = 2000000
        p = randominside(bd1)
        m1 = meancollisiontime(p, bd1, N)
        m2 = meancollisiontime(p, bd2, N)
        @test m1 ≈ m2 rtol=1e-1
    end

end

if printinfo
    println("Results:")
    println("+ Ellipse obstacle has proper arclengths.")
    println("+ randominside() works with Ellipse.")
    println("+ relocate(), collisiontime(), resolvecollision() work for")
    println("  standard particle and Ellipse.")
    println("+ Ray-splitting works with Ellipse.")
    println("+ Ellipse with equal axis is equivalent with Sinai.")
    println("+ Required time: $(round(time()-tim, digits=3)) sec.")
end
return
end
