using Test
using TEDOPA, QuadGK

@testset "TEDOPA with truncated Ohmic spectral density" begin
    rtol = 1e-6
    xc = 4
    s = 2
    ohmic_sdf = "x^a[1]"
    chainlength = 1_000

    dict = Dict(
        "environment" => Dict(
            "spectral_density_parameters" => [s],
            "spectral_density_function" => ohmic_sdf,
            "domain" => [0.0, xc],
            "temperature" => 0,
        ),
        "chain_length" => chainlength,
        "PolyChaos_nquad" => 2*chainlength,
    )

    env = chainmapping_tedopa(dict)
    sysenv_coupling = sqrt(first(quadgk(env, domain(env)...)))
    # (Here we also test the ability to directly call `env` as if it were a function.)
    @test_throws DomainError env(-1.0)

    # Exact values taken from Chin, Rivas, Huelga and Plenio, J. Math. Phys. 51 (9) 092109.
    @test frequencies(env) ≈
        [xc/2 * (1 + s^2 / ((s + 2n)*(s+2n+2))) for n in 0:(chainlength - 1)] rtol=rtol
    @test couplings(env) ≈ [
        sysenv_coupling;
        [
            xc*(1+n)*(1+s+n)/((s+2+2n)*(3+s+2n)) * sqrt((3+s+2n)/(1+s+2n)) for
            n in 0:(chainlength - 2)
        ]
    ] rtol=rtol
end

@testset "T-TEDOPA" begin
    atol = 2e-3
    x0 = 2
    xmax = 10*x0
    s = 2
    ohmic_sdf = "(x / a[1])^a[2] * exp( -x / a[1] )"
    chainlength = 10_000

    dict = Dict(
        "environment" => Dict(
            "spectral_density_parameters" => [x0, s],
            "spectral_density_function" => ohmic_sdf,
            "domain" => [0.0, xmax],
            "temperature" => 1,
        ),
        "chain_length" => chainlength,
        "PolyChaos_nquad" => 2*chainlength,
    )

    env = chainmapping_ttedopa(dict)

    @test abs(last(frequencies(env))) < atol
    @test last(couplings(env)) ≈ xmax / 2 atol=atol
end
