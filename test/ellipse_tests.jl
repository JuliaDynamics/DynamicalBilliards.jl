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




end#testset
if printinfo
    println("Results:")
    println("+ Ellipse obstacle has proper arclengths.")
    # println("+ randominside() works with Ellipse.")
    # println("+ relocate(), collisiontime(), resolvecollision() work for")
    # println("  standard particle and Ellipse.")
    # println("+ Ray-splitting works with Ellipse.")
    println("+ Required time: $(round(time()-tim, digits=3)) sec.")
end
return
end
