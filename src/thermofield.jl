function issingleton(domain)
    return minimum(domain) == maximum(domain)
end

"""
    chainmapping_thermofield(file::IOStream)

Return the frequency and coupling coefficients of the TEDOPA chain obtained by the
environment described in the JSON file `file`, after transforming it into a
`T=0` and `μ=0` environment through the thermofield procedure.

See [`chainmapping_tedopa`](@ref) for more information.
"""
function chainmapping_thermofield(file::IOStream)
    s = read(file, String)
    p = JSON.parse(s)
    parameters = merge(p, Dict("filename" => file))
    return chainmapping_thermofield(parameters)
end

"""
    chainmapping_thermofield(filename::AbstractString)

Return the frequency and coupling coefficients of the TEDOPA chain obtained by the
environment described in the JSON file called `filename`, after transforming it into a
`T=0` and `μ=0` environment through the thermofield procedure.

See [`chainmapping_tedopa`](@ref) for more information.
"""
function chainmapping_thermofield(filename::AbstractString)
    return open(filename, "r") do inputfile
        chainmapping_thermofield(inputfile)
    end
end

function prep_environments(environment::Dict{<:AbstractString,Any})
    tmp = eval(Meta.parse("(a, x) -> " * environment["spectral_density_function"]))
    sdfs = [x -> Base.invokelatest(tmp, environment["spectral_density_parameters"], x)]
    Ts = [environment["temperature"]]
    μs = [environment["chemical_potential"]]
    domains = [environment["domain"]]
    return sdfs, Ts, μs, domains
end

function prep_environments(environments::Vector{Dict{<:AbstractString,Any}})
    sdfs = []
    for d in environments
        tmp = eval(Meta.parse("(a, x) -> " * d["spectral_density_function"]))
        push!(
            sdfs, x -> Base.invokelatest(tmp, environment["spectral_density_parameters"], x)
        )
    end
    Ts = [d["temperature"] for d in environments]
    μs = [d["chemical_potential"] for d in environments]
    domains = [d["domain"] for d in environments]
    return sdfs, Ts, μs, domains
end

"""
    chainmapping_thermofield(parameters::Dict{<:AbstractString, Any})

Return the frequency and coupling coefficients of the TEDOPA chain obtained by the
environment specified by the `parameters` dictionary, after transforming it into a
`T=0` and `μ=0` environment through the thermofield procedure.

See [`chainmapping_tedopa`](@ref) for more information.
"""
function chainmapping_thermofield(parameters::Dict{<:AbstractString,Any})
    chain_length = parameters["chain_length"]

    # Is parameters["environment"] a Dict{String, Any} or a Vector{Dict{String, Any}}?
    # How do we check it? We don't! We let Julia do it by passing it to a function.
    sdfs, Ts, μs, domains = prep_environments(parameters["environment"])

    domains_empty = []
    domains_filled = []
    for (T, μ, domain) in zip(Ts, μs, domains)
        if T == 0
            # Even though Julia is able to handle T = 0 in the formulae, we still need to
            # intervene to manually restrict the domain so that the part where the
            # transformed spectral densities are identically zero are removed.
            # It might happen that some of the domains become singleton, e.g. when
            # the original domain is [0, 1] and the associated chemical potential is 0.
            # This is fine for now, we will check this later.
            push!(domains_empty, [μ, filter(>(μ), domain)...])
            push!(domains_filled, [filter(<(μ), domain)..., μ])
        else
            push!(domains_empty, domain)
            push!(domains_filled, domain)
        end
    end

    n(β, μ, ω) = (exp(β * (ω - μ)) + 1)^(-1)

    function merged_sdfempty(ω)
        return sum([(1 - n(1 / T, μ, ω)) * f(ω) for (f, T, μ) in zip(sdfs, Ts, μs)])
    end
    function merged_sdffilled(ω)
        return sum([n(1 / T, μ, ω) * f(ω) for (f, T, μ) in zip(sdfs, Ts, μs)])
    end

    # To merge the domains: concatenate, sort, then remove duplicates
    # TODO: check if there are gaps in the resulting merged domains.
    merge_domains(domains) = unique(sort(vcat(domains...)))
    merged_filled_domains = merge_domains(domains_filled)
    merged_empty_domains = merge_domains(domains_empty)

    # (We assume that the energies in the free environment Hamiltonians have already
    # been shifted with the chemical potentials. We won't do that here.)

    if !issingleton(merged_empty_domains) && !issingleton(merged_filled_domains)
        (freqempty, coupempty, sysintempty) = chainmapping(
            merged_sdfempty,
            merged_empty_domains,
            chain_length - 1;
            Nquad=parameters["PolyChaos_nquad"],
        )
        (freqfilled, coupfilled, sysintfilled) = chainmapping(
            merged_sdffilled,
            merged_filled_domains,
            chain_length - 1;
            Nquad=parameters["PolyChaos_nquad"],
        )

        chain_coefficients = Dict(
            :empty => (couplings=[sysintempty; coupempty], frequencies=freqempty),
            :filled => (couplings=[sysintfilled; coupfilled], frequencies=freqfilled),
        )
    elseif issingleton(merged_empty_domains)
        (freqfilled, coupfilled, sysintfilled) = chainmapping(
            merged_sdffilled,
            merged_filled_domains,
            chain_length - 1;
            Nquad=parameters["PolyChaos_nquad"],
        )

        chain_coefficients = Dict(
            :empty =>
                (couplings=zero([sysintfilled; coupfilled]), frequencies=zero(freqfilled)),
            :filled => (couplings=[sysintfilled; coupfilled], frequencies=freqfilled),
        )
    elseif issingleton(merged_filled_domains)
        (freqempty, coupempty, sysintempty) = chainmapping(
            merged_sdfempty,
            merged_empty_domains,
            chain_length - 1;
            Nquad=parameters["PolyChaos_nquad"],
        )

        chain_coefficients = Dict(
            :empty => (couplings=[sysintempty; coupempty], frequencies=freqempty),
            :filled =>
                (couplings=zero([sysintempty; coupempty]), frequencies=zero(freqempty)),
        )
    else  # Both merged domains are singletons. There is no output.
        error("Both merged domains are empty. Please check the input spectral densities.")
    end
    return chain_coefficients
end
