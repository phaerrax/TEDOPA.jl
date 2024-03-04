# TEDOPA

[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

TEDOPA is a package for computing the chain mapping of spectral densities of
both bosonic and fermionic environments.
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

## Environments description

Information about an environment (or a collections of environments) must be
encoded in a JSON dictionary with the following structure:

```json
{
    "environment": {
        "parameters": [1, 0.5],
        "spectral_density": "1/(10pi) * sqrt((2a[2]-a[1]+x)*(2a[2]+a[1]-x))",
        "domain": [0, 2],
        "temperature": 0.4,
        "chemical_potential": 0.5
    },
    "number_of_oscillators": 200,
    "PolyChaos_nquad": 5000
}
```

Each element of the `"environments"` subdictionary defines an environment as
follows:

* `"spectral_density"` is a string containing a valid Julia function of the
  variable `x` and possibly a set of parameters `a[1]`, `a[2]`, etc.;
* `"parameters"` is a list of numbers that will replace `a[1]`, `a[2]` and so
  on in the spectral density function;
* `"domain"` is the domain of the spectral density function;
* `"temperature"` is the temperature of the environment;
* `"chemical_potential"` is the chemical potential of the environment.

There are then two additional parameters:

* `"number_of_oscillators"` specifies the lengths of the resulting chains;
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
  mapping of the
* `chainmapping_thermofield_merge`: like the previous point, but merging
  multiple environments together (in this case, the `"environment"`
  subdictionary is replaced by a list of dictionaries, such as in the following
  example).

```json
{
    "environment": [
        {
            "parameters": [],
            "spectral_density": "1/(pi*200)",
            "domain": [-100, 100],
            "temperature": 0,
            "chemical_potential": 0.01
        },
        {
            "parameters": [],
            "spectral_density": "1/(pi*200)",
            "domain": [-100, 100],
            "temperature": 0,
            "chemical_potential": -0.01
        }
    ],
    "number_of_oscillators": 200,
    "PolyChaos_nquad": 5000
}
```
