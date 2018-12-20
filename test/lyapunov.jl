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

function test_lyapunov_values(args...)
    t = 5000.0
    radius = 1.0

    spaces = [2.0, 2.5, 3.0, 3.5, 4.0 ]
    # based on Gaspard et al (see docs) & DynamicalBilliards v2.5
    expected_values = [3.6, 1.4, 0.8, 0.6, 0.5]
    error_level = 0.1

    for (i, space) in enumerate(spaces)
        bd = billiard_polygon(6, space/(sqrt(3)); setting = "periodic")
        disc = Disk([0., 0.], radius)
        billiard = Billiard(bd.obstacles..., disc)
        p = randominside(billiard)
        λ = lyapunovspectrum(p, billiard, t)[1]
        @test abs(λ - expected_values[i]) < error_level
    end
end

billiards_testset("Lyapunov numerical values",
                  identity, caller=test_lyapunov_values)


using LinearAlgebra
function test_perturbationgrowth(p, bd, n)
    tmax = 100.0

    # there should only be floating point precision errors as we're computing
    # exactly the same thing twice
    error_level = 1e-10

    for i ∈ 1:n
        t, R, o = perturbationgrowth(p, bd, tmax)
        λ = lyapunovspectrum(p, bd, tmax)
        Δ = pertubationevolution(R)

        λ_estimate = log(norm(Δ[end]))/t[end]
        @test abs(λ[1] - λ_estimate) < error_level

        # automatically generate new intial conditions! whoopee!
        evolve!(p, bd, tmax)
    end
end

billiards_testset("Compare Lyapunovs and perturbation growth",
                  test_perturbationgrowth, 10; caller=ergodic_tests)
