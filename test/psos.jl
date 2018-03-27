using DynamicalBilliards
using Base.Test

function stadium_psos(partnum=10; printinfo = true)
    tim = time()
    @testset "PSOS: Stadium Billiard" begin
        t = 500
        l = w = 1.0
        bt = billiard_bunimovich(l, w)
        ξs, φs = poincaresection(bt, t, 5*partnum)
        #using 5*partnum so it doesn't fail all the time

        p = (10,10)
        A = zeros(Bool, p)
        c = 0
        ε = ( (2*l + w*π)/p[1], π/p[2] )
        for point ∈ zip(ξs, φs)
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
            + billiard_bunimovich is fills its PSOS
            + Required time: $(round(time()-tim, 3)) sec
            """)
        end
    end
end
