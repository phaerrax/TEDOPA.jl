# TEDOPA

[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://phaerrax.github.io/TEDOPA.jl/dev/)

TEDOPA is a package for computing the chain mapping of spectral densities of
both bosonic and fermionic environments.
Starting from the spectral density function, the temperature and the chemical
potential (if applicable), it computes the chain coefficients of the
TEDOPA-transformed environment, applying a thermalisation procedure such as
T-TEDOPA or the thermofield transformation.

## Installation

### From a registry

This package is registered in the
[TensorNetworkSimulations](https://github.com/phaerrax/TensorNetworkSimulations)
registry. If you haven't already done so, add it to your Julia installation by
running

```julia
using Pkg
pkg"registry add https://github.com/phaerrax/TensorNetworkSimulations.git"
```

(this must be done just once per Julia installation). The package can then be
installed as a normal one:

```julia
using Pkg
pkg"add TEDOPA"
```

### From GitHub

Alternatively, straight installation from GitHub is also possible:

```julia
using Pkg
pkg"add https://github.com/phaerrax/TEDOPA.jl"
```
