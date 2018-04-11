using DynamicalBilliards
using Base.Test

function stadium_psos(partnum=10; printinfo = true)
    tim = time()
    @testset "PSOS: Stadium Billiard" begin
        t = 500
        l = w = 1.0
        bt = billiard_bunimovich(l, w)
        ξs, φs = boundarymap(bt, t, 5*partnum)
        #using 5*partnum so it doesn't fail all the time

        p = (10,10)
        A = zeros(Bool, p)
        c = 0
        ε = ( (2*l + w*π)/p[1], π/p[2] )
        for point ∈ zip(vcat(ξs...), vcat(φs...))
            id = clamp.(ceil.(Int, (point .- (0, -π/2))./ε), (1,1), p)
            if !A[id...]
                A[id...] = true
                c += 1
                if c == p[1]*p[2]
                    break
                end
            end
        end
        @test *(A...)
        if printinfo
            println("""
            Results:
            + billiard_bunimovich fills its PSOS
            + Required time: $(round(time()-tim, 3)) sec
            """)
        end
    end
end

function cut_psos(partnum=10; printinfo = true)
    tim = time()
    @testset "Cut PSOS: Sinai" begin
        t = 1000
        bt = billiard_sinai()
        plane = InfiniteWall([0.5, 0.0], [0.5, 1.0], [1.0, 0.0])
        @testset "pinned" begin
            p = MagneticParticle(0.2, 0.5, -π/2, 1/0.3)
            a, b = psoscut(p,bt,plane, t)
            @test length(a) == length(b) == 1

            p = MagneticParticle(0.1, 0.5, -π/2, 1/0.05)
            a, b = psoscut(p,bt,plane, t)
            @test length(a) == length(b) == 0
        end

        @testset "psoscut ω = $ω" for ω ∈ [0, 0.5, 1.0]
            for i in 1:partnum
                p = ω == 0 ? randominside(bt) : randominside(bt, ω)
                a, b = psoscut(p,bt,plane, t)
                for j in 1:length(a)
                    @test a[j][1] ≈ 0.5
                    @test 0 < a[j][2] < 1
                    @test -1 < b[j][1] < 1
                    @test -1 < b[j][2] < 1
                end
            end
        end
    end
    if printinfo
        println("""
        Results:
        + Poincare section trhough cut works
        + Pinned particles correctly detected
        + positions and velocities are within correct bounds
        + Required time: $(round(time()-tim, 3)) sec
        """)
    end
end
