using Test
using DynamicalBilliards
using DynamicalBilliards.Testing

function coordinate_invariance(p, bd, args...)
    @testset "$(tag(p, bd)) Coordinate change" begin
        for t in 1:2e4
            i, tmin = bounce!(p, bd)
            ξ, sφ = to_bcoords(p, bd[i])
            pos, vel = from_bcoords(ξ, sφ, bd[i])
            @test *(isapprox.(p.pos, pos, atol=1e-4)...)
            @test *(isapprox.(p.vel, vel, atol=1e-4)...)
        end
    end
end

billiards_testset("Coordinate change invariance", coordinate_invariance; caller = omni_tests)



function cut_psos(args...)
    t = 1000
    bt = billiard_sinai(0.25)
    plane = InfiniteWall([0.5, 0.0], [0.5, 1.0], [-1.0, 0.0])
    @testset "PSOS: pinned particle" begin

        # case of pinned that crosses plane but not walls
        p = MagneticParticle(0.2, 0.5, -π/2, 1/0.3)
        a, b = psos(bt, plane, t, p)
        @test typeof(a) <: Vector{<:SVector}
        @test length(a) == length(b) == 1

        # case of pinned that doesn't cross anything
        p = MagneticParticle(0.1, 0.5, -π/2, 1/0.05)
        a, b = psos(bt,plane, t, p)
        @test typeof(a) <: Vector{<:SVector}
        @test length(a) == length(b) == 0

        # case of pinned but does cross periodic walls *and* section
        bd = billiard_sinai(;setting = "periodic")
        p = MagneticParticle(0.25, 0.25, -π/2 + π/4, 0.44*2)
        a, b = psos(bd, plane, t, p)
        @test typeof(a) <: Vector{<:SVector}
        @test length(a) == length(b) == 1
    end

    for ω ∈ [0.0, 1.0]
        p = ω == 0 ? randominside(bt) : randominside(bt, ω)
        @testset "$(tag(p, bt)) psos" begin
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
        a, b = psos(bt, plane, 1000, 10)
        for i in 1:length(a)
            @test length(a[i]) == length(b[i])
        end
    end
end

billiards_testset("PSOS", identity; caller = cut_psos)




# function fills_boundarymap(p, bd)
#     @testset "$(tag(p, bd)) Fills boundary map" begin
#         ξs, sφs = boundarymap(bd, 10000, p)
#
#         partition = (10,10) # partition size
#         A = falses(partition)
#         l = totallength(bd)
#         ε = (l/partition[1], 2/partition[2]) # box size
#         c = 0
#         for point ∈ zip(ξs, sφs)
#             id = clamp.(ceil.(Int, (point .- (0, -1))./ε), (1,1), partition)
#             if !A[id...]
#                 A[id...] = true
#                 c += 1
#                 if c == p[1]*p[2]
#                     break
#                 end
#             end
#         end
#         for a in A
#             @test a
#         end
#     end
# end
#
# billiards_testset("Fills boundary map", fills_boundarymap; caller = ergodic_tests)
#
#
#
#

#
# function boundarymap_portion_test(partnum = 10; printinfo = true)
#     tim = time()
#     @testset "Bunimovich" begin
#         t = 100000.0
#         bt = billiard_bunimovich()
#         @testset "ω = $ω" for ω in [0.0, 0.1]
#             for i in 1:min(partnum, 20)
#                 p = ω == 0 ? randominside(bt) : randominside(bt, ω)
#                 φ = π/4 * rand() # so that we never find bouncing walls
#                 p.vel = (cos(φ), sin(φ))
#                 ratio, dic = boundarymap_portion(bt,t, randominside(bt), 0.1)
#                 @test ratio == 1.0
#             end
#         end
#     end
#     @testset "Mushroom" begin
#         t = 100000.0
#         l = 1.0; r = 1.0
#         @testset "w = $w" for w ∈ [0.2, 0.4]
#
#             bt = billiard_mushroom(l, w, r)
#             @testset "regular" begin
#                 for i in 1:min(partnum, 10)
#                     p = MushroomTools.randomregular(l, w, r)
#                     ratio, dic = boundarymap_portion(bt,t, p, 0.1)
#                     trueratio =  MushroomTools.g_r_2D(l,w,r)
#                     # Only one regular particle covers very small amount of space:
#                     @test ratio < trueratio
#                 end
#             end
#             @testset "chaotic" begin
#                 for i in 1:min(partnum, 10)
#                     p = MushroomTools.randomchaotic(l, w, r)
#                     ratio, dic = boundarymap_portion(bt, t, p, 0.1)
#                     trueratio =  MushroomTools.g_c_2D(l,w,r)
#                     @test trueratio - 0.1 ≤ ratio ≤ trueratio + 0.1
#                 end
#             end
#
#         end
#     end
#     if printinfo
#         println("""
# Results:
# + `boundarymap_portion` works
# + Mushroom boundary map ratios are replicated correctly
# + Bunimovich stadium always gives ratio of 1.0
# + Required time: $(round(time()-tim, digits=3)) sec
# """)
#     end
# end
#
#
# function phasespace_portion_test(partnum = 10; printinfo = true)
#     tim = time()
#     @testset "Bunimovich" begin
#         t = 100000.0
#         bt = billiard_bunimovich()
#         @testset "ω = $ω" for ω in [0.0, 0.1]
#             for i in 1:min(partnum, 20)
#                 p = ω == 0 ? randominside(bt) : randominside(bt, ω)
#                 φ = π/4 * rand() # so that we never find bouncing walls
#                 p.vel = (cos(φ), sin(φ))
#                 ratio = phasespace_portion(bt,t, randominside(bt), 0.1)
#                 @test ratio == 1.0
#             end
#         end
#     end
#     @testset "Mushroom" begin
#         t = 100000.0
#         l = 1.0; r = 1.0
#         @testset "w = $w" for w ∈ [0.2, 0.4]
#             bt = billiard_mushroom(l, w, r)
#             @testset "regular" begin
#                 for i in 1:min(partnum, 10)
#                     p = MushroomTools.randomregular(l, w, r)
#                     ratio = phasespace_portion(bt,t, p, 0.1)
#                     trueratio =  MushroomTools.g_r_3D(l,w,r)
#                     # Only one regular particle covers very small amount of space:
#                     @test ratio < trueratio
#                 end
#             end
#             @testset "chaotic" begin
#                 for i in 1:min(partnum, 10)
#                     p = MushroomTools.randomchaotic(l, w, r)
#                     ratio = phasespace_portion(bt, t, p, 0.1)
#                     trueratio =  MushroomTools.g_c_3D(l,w,r)
#                     @test trueratio - 0.1 ≤ ratio ≤ trueratio + 0.1
#                 end
#             end
#
#         end
#     end
#     if printinfo
#         println("""
# Results:
# + `phasespace_portion` works
# + Mushroom phase space ratios are replicated correctly
# + Bunimovich stadium always gives ratio of 1.0
# + Required time: $(round(time()-tim, digits=3)) sec
# """)
#     end
# end
