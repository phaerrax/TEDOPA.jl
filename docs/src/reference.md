# Reference

## How to describe environments

Each variant of TEDOPA starts from the description of an environment: in this
package, information about an environment must be encoded in a dictionary
with the following structure: here is a breakdown of all the required keys.

The `"environment"` key collects the information about the environment, and
contains the following sub-keys.

* `"spectral_density_function"` is a string containing a valid Julia function
  of the variable `x`, and possibly a set of additional parameters that must be
  named `a[1]`, `a[2]`, and so on.
* `"spectral_density_parameters"` is a list of numbers that will replace `a[1]`,
  `a[2]` etc., in the given order, in the spectral density function. If no
  additional parameters are needed in the spectral density function, assign
  the empty list `[]` to this key.
* `"domain"` is the support of the spectral density function. For the purposes of
  this package, this is always an interval of the real line, than can be given
  as a list of two real numbers, representing the endpoints of the interval. If
  the spectral density function contains singular points, it is strongly advised
  to add those points to the domain (following the proper numerical order).
  These intermediate points will be taken into account by the numerical
  integration routines, leading to more precise results.  For example, for a
  spectral density function defined on ``[0,1]`` with singular points in ``1/4``
  and ``1/2`` we should set the `domain` key to `[0,0.25,0.5,1]`.
* `"temperature"` is the temperature of the environment.
* `"chemical_potential"` is the chemical potential of the environment, which is
  a number that must lie within the domain of the spectral density function.

The last two items can be omitted if the algorithm doesn't need them.

There are then two additional parameters, which are not related to the original
environments.

* `"chain_length"` specifies the lengths of the resulting chains, i.e. how many
  sites they will have.
* `"PolyChaos_nquad"`, an integer, is a parameter that will be passed to the
  numerical integration routines used to calculate the coefficients (usually a
  higher number results in a higher precision).

### Examples

The following dictionary describes a (simplified) super-Ohmic spectral density
on ``[0,2]`` with a hard cutoff: ``J(x) = x^2``. This can be encoded in a
dictionary in the aforementioned format as follows:

```julia
Dict(                                 
    "environment" => Dict(
        "spectral_density_parameters" => [],
        "spectral_density_function" => "x^2",
        "domain" => [0, 2],    
    ),
    "chain_length" => 100,
    "PolyChaos_nquad" => 200,
)
```

This other dictionary describes an environment with a semi-elliptical spectral
density

```math
J(x) = \frac{1}{2\pi} \sqrt{x(2-x)}
```

on the domain ``[0,2]``, and given temperature and chemical potential:

```julia
Dict(
    "environment" => {
        "spectral_density_parameters" => [1, 0.5],
        "spectral_density_function" => "1/2pi * sqrt((2a[2]-a[1]+x)*(2a[2]+a[1]-x))",
        "domain" => [0, 2],
        "temperature" => 0.4,
        "chemical_potential" => 0.5
    },
    "chain_length" => 200,
    "PolyChaos_nquad" => 5000,
)
```

Note that the use of spectral function parameters is not mandatory, and is
purely for convenience.

## Chain-mapping functions

This package offers several ways of computing the chain mapping, all derived
from the original TEDOPA algorithm.

* `chainmapping_tedopa`: standard chain mapping of a single environment;
* `chainmapping_ttedopa`: single-chain thermalised chain mapping for bosonic
  environments;
* `chainmapping_tftedopa`: single-chain thermalised chain mapping for fermionic
  environments;
* `chainmapping_thermofield`: a thermofield transformation and then a chain
  mapping of the resulting (fermionic) environment; it can also merge multiple
  multiple environments together (in this case, the `"environment"`
  subdictionary is replaced by a list of dictionaries, as illustrated in the
  example below).

The dictionary can be directly supplied as an argument, directly
as a Julia `Dict` object (more precisely of the `Dict{<:AbstractString,Any}`
type) to these functions.
Alternatively, the dictionary can be stored in a file, in JSON format; then you
can supply to these functions either a string containing the file name, or
the file stream object directly.

### Example

Given a Julia dictionary `dict` such as the two ones in the [How to describe
environments](@ref), we can call one of the four methods above with `dict` as an
argument, for example `chainmapping_ttedopa(dict)`.
Say we want, instead, to store the information in a JSON file called
`ohmic.json`, with the following contents:

```json
{
    "environment": {
        "spectral_density_parameters": [1, 0.5],
        "spectral_density_function": "1/2pi * sqrt((2a[2]-a[1]+x)*(2a[2]+a[1]-x))",
        "domain": [0, 2],
        "temperature": 0.4,
        "chemical_potential": 0.5
    },
    "chain_length": 200,
    "PolyChaos_nquad": 5000
}
```

We can then directly call `chainmapping_ttedopa("ohmic.json")`, or first open
the file in Julia and then pass the stream object to the chain-mapping function:

```julia
open("ohmic.json", "r") do f
    chainmapping_ttedopa(f)
end
```

### Core algorithm

This is the chain-mapping algorithm, detailed in [Chin2010](@cite), that lies at
the core of the package.
 
Let ``f \colon X \to [0, +\infty)`` be a measurable function, and
``\{p_n\}_{n=0}^{+\infty}`` be the set of monic orthogonal polynomials with
respect to the measure ``\dd f(x) = f(x)\,\dd x``. They satisfy

```math
x p_n(x) = p_{n+1}(x) + \alpha_n p_n(x) + \beta_n p_{n-1}(x)
```

for all ``n \in \N``, with ``p_{-1} \equiv 0``.

The (internal) `chainmapping` method, given the function ``f``, returns the
sequences of the _recurrence coefficients_ ``\alpha_n`` and ``\beta_n``.
Remember that, for the purposes of this package, we always assume that the
support of ``f`` is an interval on the real line.

## Standard TEDOPA chain mapping

The chain-mapping algorithm is used by TEDOPA in order to transform the original
continuous environment into a discrete chain of modes.
Consider an environment described by the free Hamiltonian

```math
H\env = \int_X x \adj{f_x} f_x \,\dd x
```

where ``f_x`` and ``\adj f_x`` are the annihilation and creation operators,
respectively, of the environment mode at frequency/energy ``x``, and, for
example, the interaction Hamiltonian

```math
H\inter = \int_X (\adj{A\sys} f_x + \adj{f_x} A\sys) \sqrt{J(x)} \,\dd x
```

although other interaction forms can be used, as long as the operator is linear
in ``f_x`` and ``\adj f_x``. 
Then, the unitary transformation

```math
b_n = \int_X \conj{u_n(x)} f_x \,\dd x, \qquad u_n(x) = \sqrt{J(x)}
\frac{p_n(x)}{\norm{p_n}}
```

generates a new set of canonical operators ``b_n`` that satisfy the same
statistics (bosonic or fermionic) as the original ``f_x`` operators. In terms of
the new mode operators,

```math
H\env = \sum_{n=0}^{+\infty} \bigl[\omega_n \adj{b_n} b_n + \kappa_n (\adj{b_n}
b_{n+1} + \adj{b_{n+1}} b_n) \bigr],
```

and

```math
H\inter = \eta (\adj{A\sys} b_0 + \adj{b_0} A\sys)
```

where ``\omega_n = \alpha_n``, ``\kappa_n = \sqrt{\beta_{n+1}}`` for ``n \in
\N`` and ``\eta = \norm{p_0}``, where the monic orthogonal polynomials and the
recurrence coefficients are taken with respect to the measure induced by ``J``.
The `chainmapping_tedopa` method takes the spectral density function ``J`` and
computes the _frequencies_ ``\omega_n``, the _couplings_ ``\kappa_n`` and the
system-environment coupling ``\eta``.

```@docs
chainmapping_tedopa
```

The function returns a struct `env` that contains the following data.

* A list of the chain mode frequencies ``\omega_n``, that can be accessed with
  `frequencies(env)`.
* A list of the chain mode couplings, whose first element is ``\eta`` followed
  by the ``\kappa_n`` coupling coefficients, and that can be retrieved with
  `couplings(env)`.
* The domain of the spectral density function, `domain(env)`.
* The spectral density function, that can be accessed simply by calling `env`
  itself with a real number as an argument.

The last two features are clearly useless in this case, as we already know the
spectral density function from the beginning, but they can be useful with the
other TEDOPA methods that transform the function.

## Thermalised TEDOPA

When the environment starts from a thermal state at inverse temperature
``\beta`` and the interaction is of the form

```math
H\inter = \int_X A\sys (f_x + \adj{f_x}) \sqrt{J(x)} \,\dd x
```

then the _thermalised TEDOPA_ approach may be used, in which we replace the
original environment with an equivalent one with the following features.

* The new environment starts from the vacuum state.
* Assuming the original domain is ``[0,x\Max]``, the new environment is
  described by a continuum of modes ``g_x`` (following the same statistics as
  the ``f_x`` modes) with frequency/energy ``x \in [-x\Max, x\Max]`` and
  Hamiltonian operators

```math
H\env' = \int_{-x\Max}^{x\Max} x \adj{g_x} g_x \,\dd x
```

and

```math
H\inter' = \int_{-x\Max}^{x\Max} A\sys (g_x + \adj{g_x}) \sqrt{J'(x)} \,\dd x,
```

where ``J'`` is the _thermalised_ spectral density function

```math
J'(x) = \frac12 \sgn x\, J(\abs{x}) \Bigl( 1 + \coth\frac{\beta x}{2} \Bigr)
```

for a bosonic environment, or

```math
J'(x) = \frac12 J(\abs{x}) \Bigl( 1 + \tanh\frac{\beta x}{2} \Bigr)
```

for a fermionic one.

The `chainmapping_ttedopa` and `chainmapping_tftedopa` methods, for bosonic and
fermionic environments respectively, return the chain coefficients and the
modified spectral density function of the thermalised environment.

```@docs
chainmapping_ttedopa
chainmapping_tftedopa
```

## Thermofield+TEDOPA (for fermionic environments)

In the more common case of fermionic environment with an exchange-type
interaction, the above thermalisation procedure is not applicable. There is
another approach, based on the _thermofield_ transformation: we introduce a
second, ancillary environment, then we shift the dependence of the initial state
of the environment on the temperature and the chemical potential to the
coupling between the system and the two environments. Namely, we transform

```math
H\env = \int_X (x-\mu) \adj{f_x} f_x \,\dd x,\\
H\inter = \int_X (\adj{A\sys} f_x + \adj{f_x} A\sys) \sqrt{J(x)} \,\dd x
```

into

```math
H\env' = \int_X x (\adj{f_{0,x}} f_{0,x} + \adj{f_{1,x}} f_{1,x}) \,\dd x,\\
H\inter' = \int_X (\adj{A\sys} f_{0,x} + \adj{f_{0,x}} A\sys) \sqrt{J_0(x)} \,\dd x + \int_X (\adj{A\sys} f_{1,x} + \adj{f_{1,x}} A\sys) \sqrt{J_1(x)} \,\dd x
```

with the two new spectral density functions

```math
J_0(x) = \biggl( 1 + \frac{1}{e^{\beta(x-\mu)} + 1} \biggr) J(x),\\
J_1(x) = \frac{1}{e^{\beta(x-\mu)} + 1} J(x).
```

The ``f_{0,x}`` modes will start from the empty state, while the ``f_{1,x}``
ones will start from the filled state.  The `chainmapping_thermofield` method
returns a dictionary with two keys, `:empty` and `:filled`, each one pointing to
the chain coefficients and the modified spectral density function of the
thermalised environment starting from the empty or filled state, respectively.

### Multiple initial environments

The thermofield+TEDOPA method can be also be applied when the system interacts
with multiple environments at once.
In order to use the algorithm in this situation, simply replace the
`"environment"` dictionary with a _list of dictionaries_, one for each
individual environment, as in the following example.

```json
{
    "environment": [
        {
            "spectral_density_parameters": [],
            "spectral_density_function": "1/(pi*200)",
            "domain": [-100, 100],
            "temperature": 0,
            "chemical_potential": 0.01
        },
        {
            "spectral_density_parameters": [],
            "spectral_density_function": "1/(pi*200)",
            "domain": [-100, 100],
            "temperature": 0,
            "chemical_potential": -0.01
        }
    ],
    "chain_length": 200,
    "PolyChaos_nquad": 5000
}
```

```@docs
chainmapping_thermofield
```
