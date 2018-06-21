using DynamicalBilliards
using Test

function coordinates_test(partnum = 500; printinfo = true)
    tim = time()
    @testset "Coordinate changes" begin
        @testset "Sinai" begin
            bd = billiard_sinai()
            aints = arcintervals(bd)
            for i in 1:partnum
                p = randominside(bd)
                i, tmin = bounce!(p, bd)
                ξ, sφ = to_bcoords(p, bd[i])
                pos, vel = from_bcoords(ξ, sφ, bd[i])

                #print differences for debug reasons
                println("pdiff\t", p.pos - pos, "\n")
                println("vdiff\t", p.vel - vel, "\n\n")

                @test *((p.pos .≈ pos)...)
                @test *((p.vel .≈ vel)...)
            end
        end
        @testset "Stadium" begin

            bd = billiard_stadium()
            aints = arcintervals(bd)
            for i in 1:partnum
                p = randominside(bd)
                i, tmin = bounce!(p, bd)
                ξ, sφ = to_bcoords(p, bd[i])
                pos, vel = from_bcoords(ξ, sφ, bd[i])

                #print differences for debug reasons
                println("pdiff\t", p.pos - pos, "\n")
                println("vdiff\t", p.vel - vel, "\n\n")

                @test *((p.pos .≈ pos)...)
                @test *((p.vel .≈ vel)...)
            end
        end
    end

    if printinfo
        println(
            """
Results:
+ from_bcoords is the inverse function of to_bcoords
  on Sinai and stadium billiards
+ Required time: $(round(time()-tim, digits=3)) sec
""")
    end
end


function stadium_bm(partnum=10; printinfo = true)
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
    end
    if printinfo
        println("""
Results:
+ boundarymap works
+ billiard_bunimovich uniformly fills its boundarymap
+ Required time: $(round(time()-tim, digits=3)) sec
""")
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
            a, b = psos(bt, plane, t, p)
            @test length(a) == length(b) == 1

            p = MagneticParticle(0.1, 0.5, -π/2, 1/0.05)
            a, b = psos(bt,plane, t, p)
            @test length(a) == length(b) == 0
        end

        @testset "psos ω = $ω" for ω ∈ [0, 0.5, 1.0]
            for i in 1:partnum
                p = ω == 0 ? randominside(bt) : randominside(bt, ω)
                a, b = psos(bt, plane, t, p)
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
+ Required time: $(round(time()-tim, digits=3)) sec
""")
    end
end


function boundarymap_portion_test(partnum = 10; printinfo = true)
    tim = time()
    @testset "Bunimovich" begin
        t = 1000000.0
        bt = billiard_bunimovich()
        @testset "ω = $ω" for ω in [0.0, 0.1]
            for i in 1:min(partnum, 20)
                p = ω == 0 ? randominside(bt) : randominside(bt, ω)
                φ = π/4 * rand() # so that we never find bouncing walls
                p.vel = (cos(φ), sin(φ))
                ratio, dic = boundarymap_portion(bt,t, randominside(bt), 0.1)
                @test ratio == 1.0
            end
        end
    end
    @testset "Mushroom" begin
        t = 1000000.0
        l = 1.0; r = 1.0
        @testset "w = $w" for w ∈ [0.2, 0.4]

            bt = billiard_mushroom(l, w, r)
            @testset "regular" begin
                for i in 1:min(partnum, 10)
                    p = MushroomTools.randomregular(l, w, r)
                    ratio, dic = boundarymap_portion(bt,t, p, 0.1)
                    trueratio =  MushroomTools.g_r_2D(l,w,r)
                    # Only one regular particle covers very small amount of space:
                    @test ratio < trueratio
                end
            end
            @testset "chaotic" begin
                for i in 1:min(partnum, 10)
                    p = MushroomTools.randomchaotic(l, w, r)
                    ratio, dic = boundarymap_portion(bt, t, p, 0.1)
                    trueratio =  MushroomTools.g_c_2D(l,w,r)
                    @test trueratio - 0.1 ≤ ratio ≤ trueratio + 0.1
                end
            end

        end
    end
    if printinfo
        println("""
Results:
+ `boundarymap_portion` works
+ Mushroom boundary map ratios are replicated correctly
+ Bunimovich stadium always gives ratio of 1.0
+ Required time: $(round(time()-tim, digits=3)) sec
""")
    end
end


function phasespace_portion_test(partnum = 10; printinfo = true)
    tim = time()
    @testset "Bunimovich" begin
        t = 1000000.0
        bt = billiard_bunimovich()
        @testset "ω = $ω" for ω in [0.0, 0.1]
            for i in 1:min(partnum, 20)
                p = ω == 0 ? randominside(bt) : randominside(bt, ω)
                φ = π/4 * rand() # so that we never find bouncing walls
                p.vel = (cos(φ), sin(φ))
                ratio = phasespace_portion(bt,t, randominside(bt), 0.1)
                @test ratio == 1.0
            end
        end
    end
    @testset "Mushroom" begin
        t = 1000000.0
        l = 1.0; r = 1.0
        @testset "w = $w" for w ∈ [0.2, 0.4]
            bt = billiard_mushroom(l, w, r)
            @testset "regular" begin
                for i in 1:min(partnum, 10)
                    p = MushroomTools.randomregular(l, w, r)
                    ratio = phasespace_portion(bt,t, p, 0.1)
                    trueratio =  MushroomTools.g_r_3D(l,w,r)
                    # Only one regular particle covers very small amount of space:
                    @test ratio < trueratio
                end
            end
            @testset "chaotic" begin
                for i in 1:min(partnum, 10)
                    p = MushroomTools.randomchaotic(l, w, r)
                    ratio = phasespace_portion(bt, t, p, 0.1)
                    trueratio =  MushroomTools.g_c_3D(l,w,r)
                    @test trueratio - 0.1 ≤ ratio ≤ trueratio + 0.1
                end
            end

        end
    end
    if printinfo
        println("""
Results:
+ `phasespace_portion` works
+ Mushroom phase space ratios are replicated correctly
+ Bunimovich stadium always gives ratio of 1.0
+ Required time: $(round(time()-tim, digits=3)) sec
""")
    end
end
