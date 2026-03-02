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

@testset "Convergence of T-TEDOPA coefficients" begin
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

@testset "Thermofield with semi-elliptical densities" begin
    rtol = 1e-6
    a = 1 / 2pi
    semielliptic_sdf = "a[1] * sqrt( x * (2-x) )"
    chainlength = 1_000

    dict_0 = Dict(
        "environment" => Dict(
            "spectral_density_parameters" => [a],
            "spectral_density_function" => semielliptic_sdf,
            "domain" => [0.0, 2.0],
        ),
        "chain_length" => chainlength,
        "PolyChaos_nquad" => 2*chainlength,
    )

    envs_0 = chainmapping_tedopa(dict_0)

    @test frequencies(envs_0) ≈ ones(chainlength) rtol=rtol
    @test couplings(envs_0) ≈ fill(0.5, chainlength) rtol=rtol

    dict_inf = Dict(
        "environment" => Dict(
            "spectral_density_parameters" => [a],
            "spectral_density_function" => semielliptic_sdf,
            "domain" => [0.0, 2.0],
            "temperature" => +Inf,
            "chemical_potential" => 1,
        ),
        "chain_length" => chainlength,
        "PolyChaos_nquad" => 2*chainlength,
    )

    envs_inf = chainmapping_thermofield(dict_inf)

    @test frequencies(envs_inf[:filled]) ≈ frequencies(envs_0) rtol=rtol
    @test couplings(envs_inf[:filled]) ≈
        [1/sqrt(2); ones(chainlength-1)] .* couplings(envs_0) rtol=rtol
    @test frequencies(envs_inf[:empty]) ≈ frequencies(envs_0) rtol=rtol
    @test couplings(envs_inf[:empty]) ≈
        [1/sqrt(2); ones(chainlength-1)] .* couplings(envs_0) rtol=rtol
end
