function issingleton(domain)
    return minimum(domain) == maximum(domain)
end

"""
    chainmapping_thermofield(parameters::Dict{<:AbstractString, Any})
    chainmapping_thermofield(file::IOStream)
    chainmapping_thermofield(filename::AbstractString)

Return the thermalised spectral densities and the associated frequency and coupling coefficients obtained by the thermofield+TEDOPA transformation of the environment specified by either
the `parameters` dictionary or by an external JSON file, passed as a stream object `file` or
by its name `filename`.
The returned object is a `NamedTuple` with two fields, `empty` and `filled`: the former contains the information relative to the final environment that starts from the vacuum state, the latter containes information about the completely filled state.

See [`chainmapping_tedopa`](@ref) for more information.

# Example

```julia-repl
julia> dict = Dict(
             "environment" => Dict(
                 "spectral_density_parameters" => [],
                 "spectral_density_function" => "1/2pi * sqrt(x*(2-x))",
                 "domain" => [0, 2],
                 "temperature" => 10.0,
                 "chemical_potential" => 0.5,
             ),
             "chain_length" => 50,
             "PolyChaos_nquad" => 200,
         );

julia> env = chainmapping_thermofield(dict);

julia> frequencies(env.empty)
50-element Vector{Float64}:
 1.0121826824869755
 1.0000050727401346
 0.999999998949
 0.9999999983153369
 0.9999999983161233
 ⋮
 0.9999999985991959
 0.9999999986176314
 0.9999999986371221
 0.9999999986577404
 0.9999999986795602

julia> couplings(env.filled)
50-element Vector{Float64}:
 0.34910972608058455
 0.49983990373980103
 0.4999998802225264
 0.4999998789916133
 0.49999984443770995
 ⋮
 0.4999984995521281
 0.4999984708186977
 0.49999844246334957
 0.49999841450784854
 0.4999983869753035

julia> domain(env.filled)
2-element Vector{Float64}:
 0.0
 2.0

julia> env.empty(0.2)
0.04703033939361449

julia> env.filled(0.2)
0.04846262646152273
```
"""
function chainmapping_thermofield end

function chainmapping_thermofield(file::IOStream)
    s = read(file, String)
    p = JSON.parse(s)
    parameters = merge(p, Dict("filename" => file))
    return chainmapping_thermofield(parameters)
end

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
    domains = [float(sort(environment["domain"]))]
    return sdfs, Ts, μs, domains
end

function prep_environments(environments::Vector{Any})
    sdfs = []
    for d in environments
        tmp = eval(Meta.parse("(a, x) -> " * d["spectral_density_function"]))
        push!(sdfs, x -> Base.invokelatest(tmp, d["spectral_density_parameters"], x))
    end
    Ts = [d["temperature"] for d in environments]
    μs = [d["chemical_potential"] for d in environments]
    domains = [float(sort(d["domain"])) for d in environments]
    return sdfs, Ts, μs, domains
end

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
        return sum([
            minimum(dom) < ω < maximum(dom) ? (1 - n(1 / T, μ, ω)) * f(ω) : zero(ω) for
            (f, T, μ, dom) in zip(sdfs, Ts, μs, domains_empty)
        ])
    end
    function merged_sdffilled(ω)
        return sum([
            minimum(dom) < ω < maximum(dom) ? n(1 / T, μ, ω) * f(ω) : zero(ω) for
            (f, T, μ, dom) in zip(sdfs, Ts, μs, domains_filled)
        ])
    end

    # To merge the domains: concatenate, sort, then remove duplicates
    # TODO: check if there are gaps in the resulting merged domains.
    merge_domains(domains) = unique(correct_minus_00, sort(vcat(domains...)))
    # We need to call `unique` with `correct_minus_00` because otherwise
    # `isequal(-0.0, 0.0)` is false, and so `unique([-0.0, 0.0]) == [-0.0, 0.0]. This leads
    # to integration errors with quadgk if the spectral density diverges in zero.

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

        chain_coefficients = (
            empty=ChainMappedEnvironment(
                merged_empty_domains, merged_sdfempty, freqempty, [sysintempty; coupempty]
            ),
            filled=ChainMappedEnvironment(
                merged_filled_domains,
                merged_sdffilled,
                freqfilled,
                [sysintfilled; coupfilled],
            ),
        )
    elseif issingleton(merged_empty_domains)
        (freqfilled, coupfilled, sysintfilled) = chainmapping(
            merged_sdffilled,
            merged_filled_domains,
            chain_length - 1;
            Nquad=parameters["PolyChaos_nquad"],
        )

        chain_coefficients = (
            empty=ChainMappedEnvironment(
                merged_empty_domains,
                zero,
                zero(freqfilled),
                zero([sysintfilled; coupfilled]),
            ),
            filled=ChainMappedEnvironment(
                merged_filled_domains,
                merged_sdffilled,
                freqfilled,
                [sysintfilled; coupfilled],
            ),
        )
    elseif issingleton(merged_filled_domains)
        (freqempty, coupempty, sysintempty) = chainmapping(
            merged_sdfempty,
            merged_empty_domains,
            chain_length - 1;
            Nquad=parameters["PolyChaos_nquad"],
        )

        chain_coefficients = (
            empty=ChainMappedEnvironment(
                merged_empty_domains, merged_sdfempty, freqempty, [sysintempty; coupempty]
            ),
            filled=ChainMappedEnvironment(
                merged_filled_domains, zero, zero(freqempty), zero([sysintempty; coupempty])
            ),
        )
    else  # Both merged domains are singletons. There is no output.
        error("Both merged domains are empty. Please check the input spectral densities.")
    end
    return chain_coefficients
end
