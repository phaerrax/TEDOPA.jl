"""
    getchaincoefficients(envparameters)

Return the Ω, η, κ coefficients of the T-TEDOPA chain, from the information given
in the `envparameters` dictionary.
"""
function getchaincoefficients(envparameters)
    n_osc = envparameters["number_of_oscillators"]
    # There are two possible cases here, according to the parameters in the JSON
    # files:
    # 1) if there is an entry "X_chain_coefficients_file" (where "X" can be "left"
    #    or "right"), then we read the coefficients from that file.
    # 2) if a specific spectral density is described through the entries
    #    "spectral_density_parameters" and "…_function", we compute the T-TEDOPA
    #    coefficients as usual.
    # A mixed-type situation in which one bath is given by a file and the other
    # by a spectral density is allowed: they are completely independent.
    # What is 𝑛𝑜𝑡 allowed, is putting in the JSON file both a coefficient file
    # and a spectral density: in this case, an error is thrown.
    if (
        haskey(envparameters, "chain_coefficients_file") &&
        !haskey(envparameters, "spectral_density_function")
    )
        ccoeffs = DataFrame(CSV.File(envparameters["chain_coefficients_file"]))
        if size(ccoeffs, 1) ≥ n_osc
            # If >, then excess coefficients are ignored.
            # If =, all coefficients are taken into account.
            # If <, there aren't enough coefficients to fill the sites, and an error
            # is thrown.
            # TODO: study the case n_osc == 1
            Ω = ccoeffs[1:n_osc, :loc]
            η = ccoeffs[1, :int]
            κ = ccoeffs[2:n_osc, :int]
        else
            error(
                "File $(envparameters["chain_coefficients_file"]) does " *
                "not contain enough coefficients.",
            )
        end
    elseif (
        !haskey(envparameters, "chain_coefficients_file") &&
        haskey(envparameters, "spectral_density_function")
    )
        # The code which creates a function from a String is taken from
        # https://stackoverflow.com/a/53134127/4160978 
        fn = envparameters["spectral_density_function"]
        tmp = eval(Meta.parse("(a, x) -> " * fn))
        sdf = x -> Base.invokelatest(tmp, envparameters["spectral_density_parameters"], x)

        ωc = envparameters["frequency_cutoff"]
        T = envparameters["temperature"]

        if T == 0
            Ω, κ, η = chainmapcoefficients(
                sdf,
                (0, ωc),
                n_osc - 1;
                Nquad=envparameters["PolyChaos_nquad"],
                discretization=lanczos,
            )
        else
            sdf_thermalised = ω -> thermalisedJ(sdf, ω, T)
            Ω, κ, η = chainmapcoefficients(
                sdf_thermalised,
                (-ωc, 0, ωc),
                n_osc - 1;
                Nquad=envparameters["PolyChaos_nquad"],
                discretization=lanczos,
            )
        end
    end
    return Ω, κ, η
end

"""
    thermalisedJ(J::Function, ω::Real, T::Real)

Create a thermalised spectral function from `J`, at the temperature
`T`, then return its value at `ω`.
"""
function thermalisedJ(J::Function, ω::Real, T::Real)
    if T == 0
        # In this case we could write f=1 and be done with it, but instead we choose
        # this more articulated way so that even if the caller doesn't already
        # exclude (-∞,0) from the support, we do it ourselves now.
        f = ω > 0 ? one(ω) : zero(ω)
    else
        f = 0.5(1 + coth(0.5 * ω / T))
    end
    return f * sign(ω) * J(abs(ω))
end

"""
    chainmapcoefficients(J::Function, support, L::Int; kwargs...)

Compute the first `L` coefficients of the chain-map Hamiltonian obtained from
the spectral density `J` defined over `support`.

# Arguments
- `J::Function`: non-negative function defined (at least) on `support`.
- `support`: support of `J`; it can be a pair of real numbers, denoting an
  interval, but can be more generally a list of increasing real numbers `[a_1, 
  a_2, a_3, …, a_N]`, that will be interpreted as ``(a_1,a_2) ∪ (a_2,a_3) ∪ ⋯ ∪ 
  (a_{N-1},a_N)``.
  The subdivision is ignored for the calculation of the chain coefficients, but
  it can become useful to compute the integral of `J` over `support`, since the
  numerical integration algorithm may need to exclude some intermediate points
  from the domain to work (e.g. if they are singularities).
- `L::Int`: number of oscillators in the chain.
- `Nquad::Int`: a parameter passed to OrthoPoly, which determines the number of
  nodes used for the quadrature method when numerical integration is performed.

Keyword arguments are passed on to the OrthoPoly constructor of the set of
orthogonal polynomials.

It returns the tuple `(Ω,κ,η)` containing
- the single-site energies `Ω`,
- the coupling coefficients `κ` between oscillators,
- the coupling coefficient `η` between the first oscillator and the system.

The spectral function `J` is associated to the recursion coefficients ``αₙ`` and
``βₙ`` which make up the recursion formula ``x πₙ(x) = πₙ₊₁(x) + αₙ πₙ(x) +
βₙ πₙ₋₁(x)`` for the monic orthogonal polynomials ``\\{πₙ\\}_{n\\in ℕ}``
determined by `J`. In the formula, ``π₋₁`` is the null polynomial.
The chain coefficients are then given by
- ``Ωᵢ = αᵢ``, with ``i\\in\\{1,\\dotsc,L\\}``,
- ``κᵢ = \\sqrt{βᵢ₊₁}``, with ``i\\in\\{1,\\dotsc,L-1\\}``,
while `η` is the integral of `J` over its support.
"""
function chainmapcoefficients(J::Function, support, L::Int; kwargs...)
    measure = PolyChaos.Measure("measure", J, (support[begin], support[end]), false, Dict())
    # We give `lanczos` as a default discretization method; if the user explicitly supplies a
    # discretization within the keyword arguments, then this default is overwritten.
    poly = PolyChaos.OrthoPoly("poly", L, measure; discretization=lanczos, kwargs...)
    # In order to build a series of L oscillators, we need the αᵢ and βᵢ
    # coefficients from i=0 to i=L-1, that means α[1:L] and β[1:L].
    # From these, the local frequencies Ωᵢ and the coupling constants κᵢ of
    # the (i,i+1) pair are given by
    #     Ωᵢ = αᵢ            for i∈{1,…,L},
    #     κᵢ = sqrt(βᵢ₊₁)    for i∈{1,…,L-1}.
    # β₀ remains unused (and rightly so, since its value is implementation-defined).
    # The oscillator-spin coupling constant η is given, lastly, by the integral
    # of J over its support.
    α = coeffs(poly)[:, 1]
    β = coeffs(poly)[:, 2]
    Ω = α
    κ = sqrt.(β[2:end])
    η = sqrt(quadgk(J, support...)[begin])
    return Ω, κ, η
end
