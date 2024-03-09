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
                    (indicator(ω, domain_neg) * J(μ - ω) + indicator(x, domain_pos) * J(μ + ω))
                else
                    0
                end
            )
    else
        newdomain = extendeddomain
        newJ =
            ω -> (
                0.5(1 + tanh(0.5ω / T)) * (
                    indicator(ω, domain_neg) * J(μ - ω) +
                    indicator(x, domain_pos) * J(μ + ω)
                )
            )
    end
    return (newdomain, newJ)
end

"""
    chainmapping_tftedopa(file::IOStream)

Return the frequency and coupling coefficients of the TEDOPA chain obtained by the
environment described in the JSON file `file`, after transforming it into a
`T=0` and `μ=0` environment through the TF-TEDOPA algorithm.

See [`chainmapping_tedopa`](@ref) for more information.
"""
function chainmapping_tftedopa(file::IOStream)
    s = read(file, String)
    p = JSON.parse(s)
    parameters = merge(p, Dict("filename" => file))
    return chainmapping_tftedopa(parameters)
end

"""
    chainmapping_tftedopa(filename::AbstractString)

Return the frequency and coupling coefficients of the TEDOPA chain obtained by the
environment described in the JSON file called `filename`, after transforming it into a
`T=0` and `μ=0` environment through the TF-TEDOPA algorithm.

See [`chainmapping_tedopa`](@ref) for more information.
"""
function chainmapping_tftedopa(filename::AbstractString)
    return open(filename, "r") do inputfile
        chainmapping_tftedopa(inputfile)
    end
end

"""
    chainmapping_tftedopa(parameters::Dict{<:AbstractString, Any})

Return the frequency and coupling coefficients of the TEDOPA chain obtained by the
environment specified by the `parameters` dictionary, after transforming it into a
`T=0` and `μ=0` environment through the TF-TEDOPA algorithm.

See [`chainmapping_tedopa`](@ref) for more information.
"""
function chainmapping_tftedopa(parameters::Dict{<:AbstractString,Any})
    chain_length = parameters["chain_length"]
    environment = parameters["environment"]

    fn = parameters["spectral_density_function"]
    tmp = eval(Meta.parse("(a, x) -> " * fn))
    sdf = x -> Base.invokelatest(tmp, environment["spectral_density_parameters"], x)

    domain = environment["domain"]
    T = environment["temperature"]
    μ = environment["chemical_potential"]
    # If [a,b] is the support of the original spectral density, the transformed one will
    # be supported on [μ-b, μ-a] ∪ [a-μ,b-μ]. This set must not contain a gap if we want the
    # new function to be in the Szegő class, so it is necessary that a ≤ μ ≤ b.
    if μ < first(domain) || μ > last(domain)
        throw(
            DomainError(
                "the spectral function: the chemical potential must lie inside the domain"
            ),
        )
    end

    cm_freqs, cm_coups, cm_syscoup = if T == 0 && μ == 0
        # Normal TEDOPA: it's simpler.
        chainmapping(sdf, domain, chain_length - 1; Nquad=parameters["PolyChaos_nquad"])
    else
        newJ, newdomain = tftedopa_sdf_transform(sdf, domain, T, μ)
        chainmapping(newJ, newdomain, chain_length - 1; Nquad=parameters["PolyChaos_nquad"])
    end

    return (frequencies=cm_freqs, couplings=[cm_syscoup; cm_coups])
end
