# [Getting Started](@id getting-started)

This guide will help you get up and running with UniDist.jl quickly.

## Installation

UniDist.jl can be installed using Julia's package manager:

```julia
using Pkg
Pkg.add("UniDist")
```

Or from the Julia REPL package mode (press `]`):

```
pkg> add UniDist
```

## Basic Usage

### Creating Distributions

Distributions are created by calling their constructor with the appropriate parameters:

```julia
using UniDist

# Normal distribution with mean=0 and standard deviation=1
normal = Normal(0.0, 1.0)

# Exponential distribution with rate=2.0
exponential = Exponential(2.0)

# Binomial distribution with n=10 trials and p=0.5 probability
binomial = Binomial(10, 0.5)

# Beta distribution with shape parameters α=2 and β=5
beta = Beta(2.0, 5.0)
```

### Computing Probabilities

#### Probability Density/Mass Function

Use `pdf` to compute the probability density (continuous) or probability mass (discrete):

```julia
# P(X = x) for discrete, f(x) for continuous
pdf(Normal(0, 1), 0.0)      # ≈ 0.3989
pdf(Binomial(10, 0.5), 5)   # ≈ 0.2461
```

#### Cumulative Distribution Function

Use `cdf` to compute P(X ≤ x):

```julia
cdf(Normal(0, 1), 1.96)     # ≈ 0.975
cdf(Exponential(1.0), 1.0)  # ≈ 0.6321
```

#### Quantile Function

Use `quantile` to find the value x such that P(X ≤ x) = p:

```julia
quantile(Normal(0, 1), 0.975)   # ≈ 1.96
quantile(Normal(0, 1), 0.5)     # = 0.0 (median)
```

### Random Sampling

Generate random samples using the `r` function:

```julia
# Generate 1000 samples from a standard normal
samples = r(Normal(0, 1), 1000)

# Generate a single sample
single = r(Normal(0, 1))
```

### Distribution Properties

```julia
d = Normal(5.0, 2.0)

mean(d)      # 5.0
var(d)       # 4.0
support(d)   # (-Inf, Inf)
```

## R-Style Shortcuts

If you're familiar with R, you can use the shorthand functions:

```julia
dist = Normal(0, 1)

d(dist, 0.0)   # Same as pdf(dist, 0.0)
p(dist, 1.96)  # Same as cdf(dist, 1.96)
q(dist, 0.975) # Same as quantile(dist, 0.975)
r(dist, 100)   # Random samples
```

## Working with Multiple Values

All functions support vectorized operations:

```julia
d = Normal(0, 1)

# Compute PDF at multiple points
pdf(d, [-2.0, -1.0, 0.0, 1.0, 2.0])

# Compute CDF at multiple points
cdf(d, [-1.96, 0.0, 1.96])

# Compute multiple quantiles
quantile(d, [0.025, 0.5, 0.975])
```

## Next Steps

- Explore the full list of [Continuous Distributions](@ref continuous-distributions) and [Discrete Distributions](@ref discrete-distributions)
- Learn about [Survival Analysis](@ref survival-analysis) functions
- See [Statistical Intervals](@ref statistical-intervals) for confidence intervals and HDI
- Check out the [Examples](@ref examples) for real-world usage scenarios
