```@meta
CurrentModule = TEDOPA
```

# TEDOPA

*This is the documentation for the TEDOPA.jl package.*

[TEDOPA](https://github.com/phaerrax/TEDOPA.jl) is a package that implements the
Time Evolving Density operator with Orthogonal Polynomials Algorithm (TEDOPA)
[Prior2010,Chin2010](@cite), a transformation that maps a continuous Gaussian
environment into a discrete chain of modes.
Starting from the spectral density function, the temperature and the chemical
potential (if applicable), the package computes the parameters that define the
TEDOPA-transformed environment.

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

## Package features

This package offers several ways of computing the chain mapping, all derived
from the original TEDOPA algorithm.

* `chainmapping_tedopa`: the standard chain mapping [Prior2010,Chin2010](@cite),
  for a single environment.
* `chainmapping_ttedopa`: thermalised chain mapping for a bosonic environment
  [Tamascelli2019](@cite).
* `chainmapping_tftedopa`: single-chain thermalised chain mapping for a
  fermionic environment [Nuesseler2020](@cite).
* `chainmapping_thermofield`: a thermofield transformation
  [deVega2015](@cite) followed by a chain mapping of the resulting (fermionic)
  environments, that can also merge multiple multiple environments together
  [Ferracin2024](@cite).

See [Reference](@ref) for a detailed explanation of the available methods.

## Bibliography

```@bibliography
```
