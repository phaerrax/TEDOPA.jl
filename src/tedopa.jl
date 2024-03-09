"""
    chainmapping_tedopa(file::IOStream)

Return the frequency and coupling coefficients of the TEDOPA chain obtained by the
environment described in the JSON file `file`.

See [`chainmapping_tedopa`](@ref) for more information.
"""
function chainmapping_tedopa(file::IOStream)
    s = read(file, String)
    p = JSON.parse(s)
    parameters = merge(p, Dict("filename" => file))
    return chainmapping_tedopa(parameters)
end

"""
    chainmapping_tedopa(filename::AbstractString)

Return the frequency and coupling coefficients of the TEDOPA chain obtained by the
environment described in the JSON file called `filename`.

See [`chainmapping_tedopa`](@ref) for more information.
"""
function chainmapping_tedopa(filename::AbstractString)
    return open(filename, "r") do inputfile
        chainmapping_tedopa(inputfile)
    end
end

"""
    chainmapping_tedopa(parameters::Dict{<:AbstractString, Any})

Return the frequency and coupling coefficients of the TEDOPA chain obtained by the
environment specified by the `parameters` dictionary.

# Example
```julia-repl
julia> env = Dict(
    "domain" => [0, 2],
    "spectral_density_parameters" => [1, 0.5],
    "spectral_density_function" => "sqrt((2a[2]-a[1]+x)*(2a[2]+a[1]-x))",
);
julia> p = Dict(
    "environment" => env,
    "chain_length" => 200,
    "PolyChaos_nquad" => 5000,
);
julia> chainmapping_tedopa(p)
(
    frequencies=[1.0000000000000004, …, 1.000000000000001],
    couplings=[0.22360679783266632, …, 0.49999999955894253],
)
```
"""
function chainmapping_tedopa(parameters::Dict{<:AbstractString,Any})
    n_osc = parameters["chain_length"]
    environment = parameters["environment"]

    # The code which creates a function from a String is taken from
    # https://stackoverflow.com/a/53134127/4160978 
    fn = environment["spectral_density_function"]
    tmp = eval(Meta.parse("(a, x) -> " * fn))
    sdf = x -> Base.invokelatest(tmp, environment["spectral_density_parameters"], x)

    cm_freqs, cm_coups, cm_syscoup = chainmapping(
        sdf,
        sort(environment["domain"]),
        n_osc - 1;
        Nquad=parameters["PolyChaos_nquad"],
        discretization=lanczos,
    )
    return (frequencies=cm_freqs, couplings=[cm_syscoup; cm_coups])
end

"""
    chainmapping(J::Function, support, L::Int; kwargs...)

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
function chainmapping(J::Function, support, L::Int; kwargs...)
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
