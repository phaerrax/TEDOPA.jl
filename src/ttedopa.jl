function boson_thermal_factor(x, T)
    if T == 0
        # In this case we could write f=1 and be done with it, but instead we choose
        # this more articulated way so that even if the caller doesn't already
        # exclude (-∞,0) from the support, we do it ourselves now.
        f = x > 0 ? one(x) : zero(x)
    else
        f = 1 / 2 * (1 + coth(0.5 * x / T))
    end
    return f
end

"""
    chainmapping_ttedopa(parameters::Dict{<:AbstractString, Any})
    chainmapping_ttedopa(file::IOStream)
    chainmapping_ttedopa(filename::AbstractString)

Return the thermalised spectral density and the associated frequency and coupling
coefficients obtained by the T-TEDOPA transformation of the environment specified by either
the `parameters` dictionary or by an external JSON file, passed as a stream object `file` or
by its name `filename`.

See [`chainmapping_tedopa`](@ref) for more information.

# Example

```julia-repl
julia> dict = Dict(
             "environment" => Dict(
                 "spectral_density_parameters" => [],
                 "spectral_density_function" => "sqrt(x) * exp(-x/0.5)",
                 "domain" => [0, 2],
                 "temperature" => 2.5,
             ),
             "chain_length" => 100,
             "PolyChaos_nquad" => 200,
         );

julia> env = chainmapping_ttedopa(dict);

julia> frequencies(env)
100-element Vector{Float64}:
  2.7937503875293512e-8
  0.2377475766338636
 -0.06307177633684016
  0.07702453114697957
 -0.07951306693183548
  ⋮
  0.006497663318373198
 -0.00645973925374582
  0.006422260704823025
 -0.006387253852889501
  0.006352720782309697

julia> couplings(env)
100-element Vector{Float64}:
 2.500086760123663
 0.000376728898182609
 1.07931125845791
 0.9389853648435562
 1.1159207167585214
 ⋮
 0.9923014215388832
 1.0076699809013268
 0.9923919032035641
 1.0075820962595274
 0.9924754152732462

julia> domain(env)
3-element Vector{Float64}:
 -2.0
 -0.0
  2.0

julia> env(1.0)
0.410505041660018
```
"""
function chainmapping_ttedopa end

function chainmapping_ttedopa(file::IOStream)
    s = read(file, String)
    p = JSON.parse(s)
    parameters = merge(p, Dict("filename" => file))
    return chainmapping_ttedopa(parameters)
end

function chainmapping_ttedopa(filename::AbstractString)
    return open(filename, "r") do inputfile
        chainmapping_ttedopa(inputfile)
    end
end

function chainmapping_ttedopa(parameters::Dict{<:AbstractString,Any})
    n_osc = parameters["chain_length"]
    environment = parameters["environment"]

    fn = environment["spectral_density_function"]
    tmp = eval(Meta.parse("(a, x) -> " * fn))
    sdf = x -> Base.invokelatest(tmp, environment["spectral_density_parameters"], x)

    domain = float(sort(environment["domain"]))
    if !iszero(first(domain))
        throw(
            DomainError(
                "the spectral function: its support must be of the form [0, m] with m>0"
            ),
        )
    end
    T = environment["temperature"]

    therm_domain = unique(correct_minus_00, [reverse(-domain); domain])
    # Reverse the domain around zero and remove duplicates: for example
    #   [0, 2, 5]  -->  [-5, -2, 0, 0, 2, 5]  -->  [-5, -2, 0, 2, 5].
    # We need to call `unique` with `correct_minus_00` because otherwise
    # `isequal(-0.0, 0.0)` is false, and so `unique([-0.0, 0.0]) == [-0.0, 0.0]. This leads
    # to integration errors with quadgk if the spectral density diverges in zero.

    therm_sdf(x) = boson_thermal_factor(x, T) * sign(x) * sdf(abs(x))

    cm_freqs, cm_coups, cm_syscoup = chainmapping(
        therm_sdf,
        therm_domain,
        n_osc - 1;
        Nquad=parameters["PolyChaos_nquad"],
        discretization=lanczos,
    )
    return ChainMappedEnvironment(therm_domain, therm_sdf, cm_freqs, [cm_syscoup; cm_coups])
end
