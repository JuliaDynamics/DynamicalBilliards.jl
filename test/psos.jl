using DynamicalBilliards
using Base.Test

function stadium_psos(partnum=10; printinfo = true)
    tim = time()
    partnum = max(partnum, 100)
    @testset "PSOS: Stadium Billiard" begin
        t = 500
        l = w = 1.0
        bt = billiard_bunimovich(l, w)
        ξs, sφs = boundarymap(bt, t, partnum)

        p = (10,10)
        A = zeros(Bool, p)
        c = 0
        ε = ( (2*l + w*π)/p[1], 2/p[2] )
        for point ∈ zip(vcat(ξs...), vcat(sφs...))
            id = clamp.(ceil.(Int, (point .- (0, -1))./ε), (1,1), p)
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
            + boundarymap works
            + billiard_bunimovich uniformly fills its boundarymap
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
        plane = InfiniteWall([0.5, 0.0], [0.5, 1.0], [-1.0, 0.0])
        @testset "pinned particle" begin
            p = MagneticParticle(0.2, 0.5, -π/2, 1/0.3)
            a, b = psos(bt,plane, t, p)
            @test length(a) == length(b) == 1

            p = MagneticParticle(0.1, 0.5, -π/2, 1/0.05)
            a, b = psos(bt,plane, t, p)
            @test length(a) == length(b) == 0
        end

        @testset "psos ω = $ω" for ω ∈ [0, 0.5, 1.0]
            for i in 1:partnum
                p = ω == 0 ? randominside(bt) : randominside(bt, ω)
                a, b = psos(bt,plane, t, p)
                for j in 1:length(a)
                    @test a[j][1] ≈ 0.5
                    @test 0 < a[j][2] < 1
                    @test -1 < b[j][1] < 1
                    @test -1 < b[j][2] < 1
                end
            end
        end
        @testset "many particles" begin
            a, b = psos(bt, plane, 1000, partnum)
            for i in 1:length(a)
                @test length(a[i]) == length(b[i])
            end
        end
    end
    if printinfo
        println("""
        Results:
        + Poincare section through cut works
        + Pinned particles correctly detected
        + positions and velocities are within correct bounds
        + Required time: $(round(time()-tim, 3)) sec
        """)
    end
end
