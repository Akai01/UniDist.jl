# UniDist.jl

**UniDist.jl** is a Julia package for working with univariate probability distributions. It provides a consistent, intuitive API for computing probability density functions, cumulative distribution functions, quantiles, random sampling, and more.

## Key Features

- **70+ distributions**: 19 discrete and 54 continuous distributions
- **Vectorized parameters**: Pass arrays as distribution parameters for batch computations
- **Familiar API**: Follows Julia's Distributions.jl conventions
- **R-style shortcuts**: Optional `d`, `p`, `q`, `r` functions for users familiar with R
- **Survival analysis**: Built-in survival function, hazard, and cumulative hazard
- **Statistical intervals**: Equal-tailed intervals and highest density intervals (HDI)

## Package Overview

```julia
using UniDist

# Create a distribution
d = Normal(0.0, 1.0)

# Compute probabilities
pdf(d, 0.0)       # 0.3989...
cdf(d, 1.96)      # 0.975...

# Generate samples
samples = r(d, 1000)

# Statistical intervals
interval(d, 0.05)  # 95% confidence interval
hdi(d, 0.95)       # 95% highest density interval
```

## Documentation Structure

- **[Getting Started](@ref getting-started)**: Installation and first steps
- **[Distributions](@ref distributions)**: Complete list of supported distributions
  - [Continuous Distributions](@ref continuous-distributions)
  - [Discrete Distributions](@ref discrete-distributions)
- **[API Reference](@ref api)**:
  - [Core Functions](@ref core-functions)
  - [Survival Analysis](@ref survival-analysis)
  - [Statistical Intervals](@ref statistical-intervals)
- **[Examples](@ref examples)**: Real-world usage scenarios
  - [Basic Usage](@ref basic-usage)
  - [Vectorized Operations](@ref vectorized-operations)
  - [Statistical Analysis](@ref statistical-analysis)
  - [Bayesian Inference](@ref bayesian-inference)

## Quick Comparison with Other Packages

| Feature | UniDist.jl | Distributions.jl |
|---------|------------|------------------|
| Core API (`pdf`, `cdf`, `quantile`) | Yes | Yes |
| R-style shortcuts (`d`, `p`, `q`, `r`) | Yes | No |
| Survival functions | Yes | Partial |
| HDI intervals | Yes | No |
| Vectorized parameters | Yes | Limited |

## References

This package is inspired by and based on the comprehensive survey of univariate distribution relationships:

> Leemis, L. M., & McQueston, J. T. (2008). **Univariate Distribution Relationships**. *The American Statistician*, 62(1), 45–53. [https://doi.org/10.1198/000313008X270448](https://doi.org/10.1198/000313008X270448)

The paper provides a detailed chart showing the relationships between 76 univariate distributions, including transformations, special cases, and limiting relationships.

## License

UniDist.jl is released under the MIT License. See [LICENSE](https://github.com/Akai01/UniDist.jl/blob/main/LICENSE) for details.
