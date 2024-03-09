```@meta
CurrentModule = TEDOPA
```

# TEDOPA

[TEDOPA](https://github.com/phaerrax/TEDOPA.jl) is a package for computing the
chain mapping of spectral densities of both bosonic and fermionic environments.
Starting from the spectral density function, the temperature and the chemical
potential (if applicable), it computes the chain coefficients of the
TEDOPA-transformed environment, applying a thermalization procedure such as
T-TEDOPA or the thermofield transformation.

## Installation

The package must be installed directly from GitHub, as it is not in any public
registry. From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```julia
add "https://github.com/phaerrax/TEDOPA.jl.git"
```

## How to describe environments

Information about an environment (or a collections of environments) must be
encoded in a JSON dictionary with the following structure:

```json
{
    "environment": {
        "spectral_density_parameters": [1, 0.5],
        "spectral_density_function": "1/(10pi) * sqrt((2a[2]-a[1]+x)*(2a[2]+a[1]-x))",
        "domain": [0, 2],
        "temperature": 0.4,
        "chemical_potential": 0.5
    },
    "chain_length": 200,
    "PolyChaos_nquad": 5000
}
```

Each element of the `"environments"` subdictionary defines an environment as
follows:

* `"spectral_density_function"` is a string containing a valid Julia function
  of the variable `x` and possibly a set of parameters `a[1]`, `a[2]`, etc.;
* `"spectral_density_parameters"` is a list of numbers that will replace `a[1]`,
  `a[2]` and so on in the spectral density function;
* `"domain"` is the domain of the spectral density function;
* `"temperature"` is the temperature of the environment;
* `"chemical_potential"` is the chemical potential of the environment.

There are then two additional parameters:

* `"chain_length"` specifies the lengths of the resulting chains;
* `"PolyChaos_nquad"` is a parameter for the numerical integration used to
  calculate the coefficients (usually a higher number means a higher precision).

## How to compute the chain coefficients

This package offers several ways of computing the chain mapping, all derived
from the original TEDOPA algorithm.

* `chainmapping_tedopa`: standard chain mapping of a single environment;
* `chainmapping_ttedopa`: thermalized chain mapping for bosonic environments;
* `chainmapping_tftedopa`: single-chain thermalized chain mapping for fermionic
  environments;
* `chainmapping_thermofield`: a thermofield transformation and then a chain
  mapping of the (fermionic) environment; it can also merge multiple
  multiple environments together (in this case, the `"environment"`
  subdictionary is replaced by a list of dictionaries, such as in the following
  example).

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

Once you have a properly formatted dictionary in a JSON file, you can call one
of the above functions by providing a file object or a string containing the
file name; you can also provide the dictionary itself, directly.

## Chain-mapping functions
### Core algorithm
```@docs
chainmapping
```

### TEDOPA variants
```@docs
chainmapping_tedopa
chainmapping_ttedopa
chainmapping_tftedopa
chainmapping_thermofield
```
