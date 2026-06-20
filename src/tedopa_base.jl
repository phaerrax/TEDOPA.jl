export chainmapping

struct ChainMappedEnvironment{F<:Function}
    domain::Vector{Float64}
    spectral_density_function::F
    frequencies::Vector{Float64}
    couplings::Vector{Float64}
end

frequencies(ce::ChainMappedEnvironment) = ce.frequencies
couplings(ce::ChainMappedEnvironment) = ce.couplings
domain(ce::ChainMappedEnvironment) = ce.domain
function (ce::ChainMappedEnvironment)(x)
    if first(domain(ce)) ≤ x ≤ last(domain(ce))
        return ce.spectral_density_function(x)
    else
        errmsg = string(
            "The spectral density function is not defined for x ∉ [",
            first(domain(ce)),
            ", ",
            last(domain(ce)),
            "].",
        )
        throw(DomainError(x, errmsg))
    end
end

"""
    chainmapping_tedopa(parameters::Dict{<:AbstractString,Any})
    chainmapping_tedopa(file::IOStream)
    chainmapping_tedopa(filename::AbstractString)

Compute the frequency and coupling coefficients of the TEDOPA chain obtained by the
environment specified by either the `parameters` dictionary or by an external JSON file,
passed as a stream object `file` or by its name `filename`.

Return a struct whose data can be retrieved via the following methods:

* `frequencies` for the single-site frquencies/energies;
* `couplings` for the coupling coefficients (the first element of the list is the coupling
  coefficient between the system and the first site of the chain);
* `domain` for the domain of the spectral density function;
* the spectral density function can be accessed by calling the struct directly, as a
  function.

# Example

```julia-repl
julia> dict = Dict(
    "environment" => Dict(
        "spectral_density_parameters" => [],
        "spectral_density_function" => "x^2",
        "domain" => [0, 2],
    ),
    "chain_length" => 100,
    "PolyChaos_nquad" => 200,
);

julia> env = chainmapping_tedopa(dict);

julia> frequencies(env)
100-element Vector{Float64}:
 1.4999999999999993
 1.1666666666666665
 1.083333333333333
 1.0499999999999998
 1.0333333333333332
 ⋮
 1.0001073883161524
 1.0001051967178594
 1.0001030715316432
 1.0001010101010104
 1.0000990099009903

julia> couplings(env)
100-element Vector{Float64}:
 1.632993161855452
 0.3872983346207419
 0.45074893585520914
 0.47245559126153397
 0.48241815132442156
 ⋮
 0.49995252761390824
 0.4999535013925803
 0.4999544455135774
 0.49995536116913225
 0.49995624949217776

julia> domain(env)
2-element Vector{Float64}:
 0.0
 2.0

julia> env(2)
4
```
"""
function chainmapping_tedopa end

function chainmapping_tedopa(file::IOStream)
    s = read(file, String)
    p = JSON.parse(s)
    parameters = merge(p, Dict("filename" => file))
    return chainmapping_tedopa(parameters)
end

function chainmapping_tedopa(filename::AbstractString)
    return open(filename, "r") do inputfile
        chainmapping_tedopa(inputfile)
    end
end

function chainmapping_tedopa(parameters::Dict{<:AbstractString,Any})
    nsites = parameters["chain_length"]
    environment = parameters["environment"]
    domain = float(sort(environment["domain"]))

    # The code which creates a function from a String is taken from
    # https://stackoverflow.com/a/53134127/4160978 
    fn = environment["spectral_density_function"]
    tmp = eval(Meta.parse("(a, x) -> " * fn))
    sdf = x -> Base.invokelatest(tmp, environment["spectral_density_parameters"], x)

    cm_freqs, cm_coups, cm_syscoup = chainmapping(
        sdf, domain, nsites; Nquad=parameters["PolyChaos_nquad"], discretization=lanczos
    )
    return ChainMappedEnvironment(domain, sdf, cm_freqs, [cm_syscoup; cm_coups])
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
- `L::Int`: number of sites in the chain.
- `Nquad::Int`: a parameter passed to OrthoPoly, which determines the number of
  nodes used for the quadrature method when numerical integration is performed.

Keyword arguments are passed on to the OrthoPoly constructor of the set of
orthogonal polynomials.

It returns the tuple `(ω, κ, η)` containing
- the single-site energies `ω` of the first `L` sites,
- the coupling coefficients `κ` between the first `L` sites,
- the coupling coefficient `η` between the first chain site and the system.

The spectral function `J` is associated to the recursion coefficients ``αₙ`` and
``βₙ`` which make up the recursion formula ``x πₙ(x) = πₙ₊₁(x) + αₙ πₙ(x) +
βₙ πₙ₋₁(x)`` for the monic orthogonal polynomials ``\\{πₙ\\}_{n\\in ℕ}``
determined by `J`. In the formula, ``π₋₁`` is the null polynomial.
The chain coefficients are then given by
- ``ωᵢ = αᵢ``, with ``i\\in\\{1,\\dotsc,L\\}``,
- ``κᵢ = \\sqrt{βᵢ₊₁}``, with ``i\\in\\{1,\\dotsc,L-1\\}``,
while `η` is the integral of `J` over its support.
"""
function chainmapping(J::Function, support, L::Int; kwargs...)
    measure = PolyChaos.Measure("measure", J, (support[begin], support[end]), false, Dict())
    # We give `lanczos` as a default discretization method; if the user explicitly supplies
    # a discretization within the keyword arguments, then this default is overwritten.
    poly = PolyChaos.OrthoPoly("poly", L, measure; discretization=lanczos, kwargs...)
    # In order to build a series of L sites, we need the αᵢ and βᵢ
    # coefficients from i=0 to i=L-1, that means α[1:L] and β[1:L].
    # From these, the local frequencies ωᵢ and the coupling constants κᵢ of
    # the (i,i+1) pair are given by
    #     ωᵢ = αᵢ            for i∈{0,…,L},
    #     κᵢ = sqrt(βᵢ₊₁)    for i∈{0,…,L-1}.
    #     η = sqrt(∫ J)
    # β₀ remains unused. Technically, it is not part of the recurrence coefficients, but
    # libraries usually define it as the integral of the measure over its support.
    # This is what PolyChaos does, but it seems like it's not as precise as calling quadgk
    # directly: if we use η = sqrt(β[1]), then the tests fail, as β[1] is not equal to J's
    # integral within the error bounds provided by quadgk.

    α = coeffs(poly)[:, 1]  # There are L+1 elements here, since it goes from α_0 to α_L
    β = coeffs(poly)[:, 2]
    ω = α[1:L]
    κ = sqrt.(β[2:L])
    η = sqrt(quadgk(J, support...)[begin])
    # η = sqrt(β[1])
    return ω, κ, η
end

@inline correct_minus_00(x) = (x==0.0 ? zero(x) : x)
# Small helper function that treats -0.0 and 0.0 as equal. Useful when you do not want
# to discriminate between the two numbers, i.e. when using `unique`.
