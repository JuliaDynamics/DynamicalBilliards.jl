using Test
using DynamicalBilliards
using DynamicalBilliards.Testing


function test_lyapunov_spectrum(p, bd, t = 1e6, error_level = 1e-3)
    # calculate spectrum
    exps = lyapunovspectrum!(p, bd, t)
    sumpair = exps[1] + exps[4]

    # these properties should be true for all billiards
    @testset "λ₁ + λ₄ ≈ 0" begin
        @test abs(sumpair) < error_level
    end

    @testset "λ₂ ≈ 0, λ₃ ≈ 0" begin
        @test abs(exps[2]) < error_level
        @test abs(exps[3]) < error_level
    end    
end

billiards_testset("Properties of Lyapunov spectrum",
                  test_lyapunov_spectrum; caller = simple_tests)

function test_lyapunov_magnetic_limit(args...)
    bd = billiard_bunimovich(0.1)

    t = 1e5
    error_level = 5e-2
    
    for j ∈ 1:10
        pmag = randominside(bd, 1e-3)
        plin = Particle(pmag.pos, pmag.vel, SVector{2,Float64}(0,0))
        
        exp_mag = lyapunovspectrum!(pmag, bd, t)
        exp_lin = lyapunovspectrum!(plin, bd, t)
        
        @test exp_mag[1] - exp_lin[1] < error_level
    end 
end

billiards_testset("Lyapunov spectrum for small B",
                  identity; caller=test_lyapunov_magnetic_limit)
