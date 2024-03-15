# TEDOPA

[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

TEDOPA is a package for computing the chain mapping of spectral densities of
both bosonic and fermionic environments.
Starting from the spectral density function, the temperature and the chemical
potential (if applicable), it computes the chain coefficients of the
TEDOPA-transformed environment, applying a thermalization procedure such as
T-TEDOPA or the thermofield transformation.

## Installation

### From a registry

This package is registered in my
[TensorNetworkSimulations](https://github.com/phaerrax/TensorNetworkSimulations)
registry. By first adding this registry, with

```julia
using Pkg
pkg"registry add https://github.com/phaerrax/TensorNetworkSimulations.git"
```

(this must be done just once per Julia installation) the package can then be
installed as a normal one:

```julia
using Pkg
Pkg.add("TEDOPA")
```

### From GitHub

Alternatively, straight installation from GitHub is also possible:

```julia
using Pkg
Pkg.add("https://github.com/phaerrax/TEDOPA.jl")
```

## Environments description

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

## Calculating the chain coefficients

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
