function indicator(x, lb, rb)
    return lb <= x <= rb ? one(x) : zero(x)
end

function indicator(x, domain)
    return indicator(x, first(domain), last(domain))
end

function tftedopa_sdf_transform(J::Function, domain, T, μ)
    domain_neg = reverse(-domain) .+ μ
    domain_pos = domain .- μ
    extendeddomain = unique(sort([domain_neg; domain_pos]))

    if T == 0  # Cut off part of domain where new_J is zero (ω < 0)
        newdomain = [0; filter(>(0), extendeddomain)]
        newJ =
            ω -> (
                if ω >= 0
                    (iszero(indicator(ω, domain_neg)) ? zero(ω) : J(μ - ω)) +
                    (iszero(indicator(ω, domain_pos)) ? zero(ω) : J(μ + ω))
                else
                    0
                end
            )
    else
        newdomain = extendeddomain
        newJ =
            ω -> (
                0.5(1 + tanh(0.5ω / T)) * (
                    (iszero(indicator(ω, domain_neg)) ? zero(ω) : J(μ - ω)) +
                    (iszero(indicator(ω, domain_pos)) ? zero(ω) : J(μ + ω))
                )
            )
    end
    return (newJ, newdomain)
end

"""
    chainmapping_tftedopa(parameters::Dict{<:AbstractString, Any})
    chainmapping_tftedopa(file::IOStream)
    chainmapping_tftedopa(filename::AbstractString)

Return the thermalised spectral density and the associated frequency and coupling
coefficients obtained by the TF-TEDOPA transformation of the environment specified by either
the `parameters` dictionary or by an external JSON file, passed as a stream object `file` or
by its name `filename`.

See [`chainmapping_tedopa`](@ref) for more information.

# Example

```julia-repl
julia> dict = Dict(
             "environment" => Dict(
                 "spectral_density_parameters" => [],
                 "spectral_density_function" => "1/pi * sqrt(x*(2-x))",
                 "domain" => [0, 2],
                 "temperature" => 10.0,
                 "chemical_potential" => 1,
             ),
             "chain_length" => 100,
             "PolyChaos_nquad" => 200,
         );

julia> env = chainmapping_tftedopa(dict);

julia> frequencies(env)
100-element Vector{Float64}:
  0.012494792328180852
  5.2026993957315084e-6
 -1.077894557975692e-9
 -1.7277946434351792e-9
 -1.7269887048876722e-9
  ⋮
  7.564160055262281e-9
  8.971713063438136e-9
  1.0968847909825774e-8
  1.4329892390751997e-8
  2.447297564427131e-8

julia> couplings(env)
100-element Vector{Float64}:
 0.707106781448028
 0.499843803889265
 0.4999998810351929
 0.49999987899169407
 0.49999984443772305
 ⋮
 0.4999988694272515
 0.4999990207121836
 0.49999920014918586
 0.4999994195321602
 0.4999997061585101

julia> domain(env)
2-element Vector{Float64}:
 -1.0
  1.0

julia> env(-0.2)
0.308760037243951
```
"""
function chainmapping_tftedopa end

function chainmapping_tftedopa(file::IOStream)
    s = read(file, String)
    p = JSON.parse(s)
    parameters = merge(p, Dict("filename" => file))
    return chainmapping_tftedopa(parameters)
end

function chainmapping_tftedopa(filename::AbstractString)
    return open(filename, "r") do inputfile
        chainmapping_tftedopa(inputfile)
    end
end

function chainmapping_tftedopa(parameters::Dict{<:AbstractString,Any})
    chain_length = parameters["chain_length"]
    environment = parameters["environment"]

    fn = environment["spectral_density_function"]
    tmp = eval(Meta.parse("(a, x) -> " * fn))
    sdf = x -> Base.invokelatest(tmp, environment["spectral_density_parameters"], x)

    domain = float(sort(environment["domain"]))
    T = environment["temperature"]
    μ = environment["chemical_potential"]
    # If [a,b] is the support of the original spectral density, the transformed one will
    # be supported on [μ-b, μ-a] ∪ [a-μ,b-μ]. This set must not contain a gap if we want the
    # new function to be in the Szegő class, so it is necessary that a ≤ μ ≤ b.
    if μ < first(domain) || μ > last(domain)
        throw(
            DomainError(
                "the spectral function", "the chemical potential must lie inside the domain"
            ),
        )
    end

    return if T == 0 && μ == 0
        # Normal TEDOPA: it's simpler.
        cm_freqs, cm_coups, cm_syscoup = chainmapping(
            sdf, domain, chain_length - 1; Nquad=parameters["PolyChaos_nquad"]
        )
        ChainMappedEnvironment(domain, sdf, cm_freqs, [cm_syscoup; cm_coups])
    else
        newJ, newdomain = tftedopa_sdf_transform(sdf, domain, T, μ)
        cm_freqs, cm_coups, cm_syscoup = chainmapping(
            newJ, newdomain, chain_length - 1; Nquad=parameters["PolyChaos_nquad"]
        )
        ChainMappedEnvironment(newdomain, newJ, cm_freqs, [cm_syscoup; cm_coups])
    end
end
