# UniDist.jl

[![CI](https://github.com/Akai01/UniDist.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/Akai01/UniDist.jl/actions/workflows/CI.yml)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://Akai01.github.io/UniDist.jl/dev)

> **Note:** This package is experimental. The API may change without notice in future versions.

A Julia package for univariate probability distributions with vectorized parameter support.

## Installation

```julia
using Pkg
Pkg.add("UniDist")
```

## Quick Start

```julia
using UniDist

# Create distributions
d = Normal(0.0, 1.0)

# Core functions
pdf(d, 0.0)       # Probability density
cdf(d, 0.0)       # Cumulative distribution
quantile(d, 0.5)  # Quantile function
mean(d)           # Mean
var(d)            # Variance

# R-style shortcuts
d(dist, x)  # pdf
p(dist, x)  # cdf
q(dist, u)  # quantile
r(dist, n)  # random samples

# Survival analysis
sf(d, x)      # Survival function
hazard(d, x)  # Hazard function
cumhaz(d, x)  # Cumulative hazard

# Statistical intervals
interval(d, 0.05)  # 95% equal-tailed interval
hdi(d, 0.95)       # 95% highest density interval
```

## Supported Distributions

**Discrete:** Benford, Bernoulli, BetaBinomial, BetaPascal, Binomial, DiscreteUniform, DiscreteWeibull, GammaPoisson, Geometric, Hypergeometric, Logarithm, NegativeHypergeometric, Pascal, Poisson, Polya, PowerSeries, Rectangular, Zeta, Zipf

**Continuous:** Arcsin, Arctangent, Beta, Cauchy, Chi, ChiSquare, Erlang, Error, Exponential, ExponentialPower, ExtremeValue, F, Gamma, GammaNormal, GeneralizedGamma, GeneralizedPareto, Gompertz, HyperbolicSecant, Hyperexponential, Hypoexponential, IDB, InverseGaussian, InvertedBeta, InvertedGamma, KolmogorovSmirnov, Laplace, LogGamma, Logistic, LogisticExponential, LogLogistic, LogNormal, Lomax, Makeham, Minimax, Muth, NoncentralBeta, NoncentralChiSquare, NoncentralF, NoncentralT, DoublyNoncentralF, DoublyNoncentralT, Normal, Pareto, Power, Rayleigh, StandardCauchy, StandardNormal, StandardPower, StandardTriangular, StandardUniform, Triangular, Uniform, VonMises, Weibull

## References

This package is inspired by the comprehensive survey of univariate distribution relationships:

> Leemis, L. M., & McQueston, J. T. (2008). **Univariate Distribution Relationships**. *The American Statistician*, 62(1), 45–53. [DOI: 10.1198/000313008X270448](https://doi.org/10.1198/000313008X270448)

## License

MIT License - see [LICENSE](LICENSE) for details.
