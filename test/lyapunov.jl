using Base.Test

function lyapunov_spectrum(partnum=500; printinfo = true)
tim = time()
partnum= min(10, partnum)
@testset "Lyapunov Spectrum (straight)" begin


    @testset "Hexagonal Periodic" begin
        l = 2.0
        r = 1.0
        sides = 6
        bt = billiard_polygon(6,l; setting = "periodic")
        disc = Disk([0., 0.], r)
        push!(bt, disc)
        tt=10000.0

        for i in 1:partnum
            p = randominside(bt)
            exps = lyapunovspectrum!(p, bt, tt)
            error_level = 1e-2
            sumpair = exps[1] + exps[4]

            @test abs(sumpair) < error_level
            @test abs(exps[2]) < error_level
            @test abs(exps[3]) < error_level
        end#particle loop
    end
    @testset "Sinai" begin
        bt = billiard_sinai()
        tt=10000.0
        for i in 1:partnum
            p = randominside(bt)
            exps = lyapunovspectrum!(p, bt, tt)
            error_level = 1e-2
            sumpair = exps[1] + exps[4]

            @test abs(sumpair) < error_level
            @test abs(exps[2]) < error_level
            @test abs(exps[3]) < error_level
        end#particle loop
    end
end
if printinfo
    println("Results:")
    println("+ lyapunovspectrum() call works on")
    println("  hexagonal lorentz and sinai billiard.")
    println("+ λ₁ + λ₄ ≈ 0.")
    println("+ λ₂ ≈ λ₃ ≈ 0.")
    println("+ Required time: $(round(time()-tim, 3)) sec.")
end
return
end
