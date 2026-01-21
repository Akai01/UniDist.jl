# [Core Functions](@id core-functions)

This page documents the core API functions for working with distributions in UniDist.jl.

## Probability Functions

### `pdf` - Probability Density/Mass Function

Computes the probability density (continuous) or probability mass (discrete) at a given point.

```julia
pdf(dist, x)
```

**Arguments:**
- `dist`: A distribution object
- `x`: Point(s) at which to evaluate

**Returns:** Probability density or mass value(s)

**Examples:**
```julia
# Continuous: probability density
d = Normal(0, 1)
pdf(d, 0.0)           # ≈ 0.3989

# Discrete: probability mass
d = Binomial(10, 0.5)
pdf(d, 5)             # ≈ 0.2461

# Multiple points
pdf(Normal(0, 1), [-1.0, 0.0, 1.0])  # [0.242, 0.399, 0.242]
```

---

### `cdf` - Cumulative Distribution Function

Computes P(X ≤ x), the probability that the random variable is less than or equal to x.

```julia
cdf(dist, x)
```

**Arguments:**
- `dist`: A distribution object
- `x`: Point(s) at which to evaluate

**Returns:** Cumulative probability value(s)

**Examples:**
```julia
d = Normal(0, 1)

cdf(d, 0.0)           # 0.5 (half of probability below mean)
cdf(d, 1.96)          # ≈ 0.975

# P(X > x) = 1 - cdf(x)
1 - cdf(d, 1.96)      # ≈ 0.025

# Multiple points
cdf(d, [-1.96, 0, 1.96])  # [0.025, 0.5, 0.975]
```

---

### `quantile` - Quantile Function (Inverse CDF)

Finds the value x such that P(X ≤ x) = p.

```julia
quantile(dist, p)
```

**Arguments:**
- `dist`: A distribution object
- `p`: Probability value(s), 0 ≤ p ≤ 1

**Returns:** Quantile value(s)

**Examples:**
```julia
d = Normal(0, 1)

quantile(d, 0.5)      # 0.0 (median)
quantile(d, 0.975)    # ≈ 1.96
quantile(d, 0.025)    # ≈ -1.96

# Multiple quantiles
quantile(d, [0.025, 0.5, 0.975])  # [-1.96, 0.0, 1.96]
```

---

## R-Style Shorthand Functions

For users familiar with R's syntax, UniDist.jl provides shorthand aliases:

| Shorthand | Full Function | R Equivalent |
|-----------|---------------|--------------|
| `d(dist, x)` | `pdf(dist, x)` | `dnorm(x)` |
| `p(dist, x)` | `cdf(dist, x)` | `pnorm(x)` |
| `q(dist, p)` | `quantile(dist, p)` | `qnorm(p)` |
| `r(dist, n)` | Random samples | `rnorm(n)` |

**Examples:**
```julia
dist = Normal(0, 1)

d(dist, 0.0)          # Same as pdf(dist, 0.0)
p(dist, 1.96)         # Same as cdf(dist, 1.96)
q(dist, 0.975)        # Same as quantile(dist, 0.975)
r(dist, 100)          # 100 random samples
```

---

## Distribution Properties

### `mean` - Expected Value

Returns the mean (expected value) of the distribution.

```julia
mean(dist)
```

**Examples:**
```julia
mean(Normal(5, 2))        # 5.0
mean(Exponential(2))      # 2.0
mean(Binomial(10, 0.3))   # 3.0
mean(Poisson(4.5))        # 4.5
```

---

### `var` - Variance

Returns the variance of the distribution.

```julia
var(dist)
```

**Examples:**
```julia
var(Normal(0, 2))         # 4.0 (σ²)
var(Exponential(2))       # 4.0 (θ²)
var(Binomial(10, 0.3))    # 2.1 (np(1-p))
var(Poisson(4.5))         # 4.5 (λ)
```

---

### `support` - Distribution Support

Returns the support (range of possible values) of the distribution.

```julia
support(dist)
```

**Returns:** A tuple `(lower, upper)` for continuous distributions, or the range for discrete distributions.

**Examples:**
```julia
support(Normal(0, 1))         # (-Inf, Inf)
support(Exponential(1))       # (0, Inf)
support(Beta(2, 3))           # (0, 1)
support(Uniform(0, 10))       # (0, 10)
support(Binomial(10, 0.5))    # Range or bounds
```

---

## Random Sampling

### `r` - Generate Random Samples

Generates random samples from the distribution.

```julia
r(dist, n=1; rng=default_rng())
```

**Arguments:**
- `dist`: A distribution object
- `n`: Number of samples (default: 1)
- `rng`: Random number generator (optional)

**Returns:** Array of random samples

**Examples:**
```julia
d = Normal(0, 1)

# Single sample
r(d)                  # e.g., 0.342

# Multiple samples
samples = r(d, 1000)
length(samples)       # 1000

# With specific RNG for reproducibility
using Random
rng = MersenneTwister(42)
r(d, 5; rng=rng)      # Reproducible samples
```

---

## Vectorized Operations

All core functions support vectorized operations over multiple input values:

```julia
d = Normal(0, 1)
xs = [-2.0, -1.0, 0.0, 1.0, 2.0]

# PDF at multiple points
pdf(d, xs)            # [0.054, 0.242, 0.399, 0.242, 0.054]

# CDF at multiple points
cdf(d, xs)            # [0.023, 0.159, 0.5, 0.841, 0.977]

# Multiple quantiles
ps = [0.1, 0.25, 0.5, 0.75, 0.9]
quantile(d, ps)       # [-1.28, -0.67, 0.0, 0.67, 1.28]
```

---

## Type Hierarchy

UniDist.jl uses a type hierarchy to distinguish between distribution types:

```julia
abstract type AbstractDistribution end
abstract type DiscreteDistribution <: AbstractDistribution end
abstract type ContinuousDistribution <: AbstractDistribution end
```

This allows for type-based dispatch and specialization:

```julia
# Check distribution type
d = Normal(0, 1)
d isa ContinuousDistribution  # true
d isa DiscreteDistribution    # false

d = Binomial(10, 0.5)
d isa DiscreteDistribution    # true
```
