function boson_thermal_factor(x, T)
    if T == 0
        # In this case we could write f=1 and be done with it, but instead we choose
        # this more articulated way so that even if the caller doesn't already
        # exclude (-∞,0) from the support, we do it ourselves now.
        f = x > 0 ? one(x) : zero(x)
    else
        f = 1/2 * (1 + coth(0.5 * x / T))
    end
    return f
end

"""
    chainmapping_ttedopa(parameters::Dict{AbstractString, Any})

Return the frequency and coupling coefficients of the TEDOPA chain obtained by the
environment specified by the `envparameters` dictionary, after a thermalization procedure.

See [`chainmapping_tedopa`](@ref) for more information.
"""
function chainmapping_ttedopa(parameters::Dict{<:AbstractString,Any})
    n_osc = parameters["number_of_oscillators"]
    environment = parameters["environment"]

    fn = environment["spectral_density_function"]
    tmp = eval(Meta.parse("(a, x) -> " * fn))
    sdf = x -> Base.invokelatest(tmp, environment["spectral_density_parameters"], x)

    domain = environment["domain"]
    if first(domain) != 0
        throw(
            DomainError(
                "the spectral function: its support must be of the form [0, m] with m>0"
            ),
        )
    end
    T = environment["temperature"]
    # Reverse the domain around zero and remove duplicates: for example
    #   [0, 2, 5]  -->  [-5, -2, 0, 0, 2, 5]  -->  [-5, -2, 0, 2, 5].
    therm_domain = unique([reverse(-domain); domain])
    therm_sdf(x) = boson_thermal_factor(x, T) * sign(x) * sdf(abs(x))

    cm_freqs, cm_coups, cm_syscoup = chainmapping(
        therm_sdf,
        therm_domain,
        n_osc - 1;
        Nquad=parameters["PolyChaos_nquad"],
        discretization=lanczos,
    )
    return (frequencies=cm_freqs, couplings=[cm_syscoup; cm_coups])
end
