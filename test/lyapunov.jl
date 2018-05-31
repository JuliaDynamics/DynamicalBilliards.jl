using Base.Test

function lyapunov_spectrum(partnum=500; printinfo = true)
    tim = time()
    partnum= min(10, partnum)
    @testset "Lyapunov Spectrum (straight)" begin


        @testset "Hexagonal Periodic" begin
            l = 2.0
            r = 1.0
            sides = 6
            bt = billiard_hexagonal_sinai(r,l; setting = "periodic")
            #disc = Disk([0., 0.], r)
            #push!(bt, disc)
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
        println("+ lyapunovspectrum!() call works on")
        println("  hexagonal lorentz and sinai billiard.")
        println("+ λ₁ + λ₄ ≈ 0.")
        println("+ λ₂ ≈ λ₃ ≈ 0.")
        println("+ Required time: $(round(time()-tim, 3)) sec.")
    end
    return
end

function lyapunov_magnetic(partnum=500; printinfo = true)
    tim = time()
    partnum= min(10, partnum)
    @testset "Lyapunov Spectrum (magnetic)" begin
        @testset "Sinai billiard, ω = 0.75" begin
            bt = billiard_sinai()
            tt=10000.0
            error_level = 1e-2
            for i in 1:partnum
                p = randominside(bt, 0.75)
                exps = lyapunovspectrum!(p, bt, tt)
                sumpair = exps[1] + exps[4]

                @test abs(sumpair) < error_level
                @test abs(exps[2]) < error_level
                @test abs(exps[3]) < error_level
            end#particle loop
        end

        @testset "Sinai billiard, ω = 2" begin
            bt = billiard_sinai(;setting="periodic")
            p = MagneticParticle(0.1, 0.5, -π/2, 2.0)
            error_level = 1e-5
            Λ = lyapunovspectrum!(p, bt, 10000.0)
            @test abs(Λ[1]) < error_level
        end

        @testset "Hexagonal Sinai, ω = 0.01" begin
            bt = billiard_hexagonal_sinai(0.5, 1.0; setting = "periodic")
            tt=10000.0
            error_level = 5e-2

            emagsum = 0.0; elinsum = 0.0
            for j ∈ 1:10
                pmag = randominside(bt, 0.01)
                plin = Particle(pmag.pos, pmag.vel, SVector{2,Float64}(0,0))

                emag = lyapunovspectrum!(pmag, bt, tt)
                elin = lyapunovspectrum!(plin, bt, tt)
                emagsum += emag[1]; elinsum += elin[1]
            end #particle loop
            difference = emagsum - elinsum
            @test abs(difference) < partnum*error_level
        end

        @testset "Stadium, ω = 1e-3" begin
            bt = billiard_bunimovich(0.1)
            tt=10000.0
            error_level = 5e-2

            emagsum = 0.0; elinsum = 0.0
            for j ∈ 1:10
                pmag = randominside(bt, 1e-3)
                plin = Particle(pmag.pos, pmag.vel, SVector{2,Float64}(0,0))

                emag = lyapunovspectrum!(pmag, bt, tt)
                elin = lyapunovspectrum!(plin, bt, tt)
                emagsum += emag[1]; elinsum += elin[1]
            end #particle loop
            difference = emagsum - elinsum
            @test abs(difference) < partnum*error_level
        end

    end
    if printinfo
        println("""Results:
                + lyapunovspectrum!() call works on
                  magnetic Sinai billiards, magnetic
                  hexagonal Sinai billiards and
                  magnetic stadium billiards
                + λ₁ + λ₄ ≈ 0 and λ₂ ≈ λ₃ ≈ 0 for
                  Sinai billiards
                + λ₁ ≈ 0 for pinned particles
                + magnetic λ similar to non-magnetic
                  λ for small ω
                + Required time: $(round(time()-tim, 3)) sec.
                """)
    end
    return
end
