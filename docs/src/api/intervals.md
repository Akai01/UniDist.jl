# [Statistical Intervals](@id statistical-intervals)

UniDist.jl provides functions for computing statistical intervals from distributions, useful for uncertainty quantification, Bayesian inference, and hypothesis testing.

## Overview

| Function | Description | Best For |
|----------|-------------|----------|
| `interval(d, α)` | Equal-tailed interval | Symmetric distributions |
| `hdi(d, mass)` | Highest Density Interval | Skewed/multimodal distributions |

---

## Equal-Tailed Interval

### `interval` - Symmetric Probability Interval

Computes an equal-tailed interval that contains (1 - α) × 100% of the probability mass.

```julia
interval(dist, alpha=0.05)
```

**Arguments:**
- `dist`: A distribution object
- `alpha`: Significance level (default: 0.05 for 95% interval)

**Returns:** Tuple `(lower, upper)` where:
- `lower = quantile(dist, α/2)`
- `upper = quantile(dist, 1 - α/2)`

**Examples:**
```julia
# 95% interval for standard normal
d = Normal(0, 1)
interval(d, 0.05)     # (-1.96, 1.96)

# 90% interval
interval(d, 0.10)     # (-1.645, 1.645)

# 99% interval
interval(d, 0.01)     # (-2.576, 2.576)
```

**Properties:**
- Equal probability in each tail: P(X < lower) = P(X > upper) = α/2
- Always symmetric in probability, not necessarily in value
- Simple and widely understood

---

## Highest Density Interval

### `hdi` - Highest Density Interval (HDI/HPD)

Computes the Highest Density Interval - the narrowest interval containing the specified probability mass.

```julia
hdi(dist, mass=0.9; grid=2048, eps=1e-6)  # Continuous
hdi(dist, mass=0.9)                        # Discrete
```

**Arguments:**
- `dist`: A distribution object
- `mass`: Probability mass to contain (default: 0.9 for 90% HDI)
- `grid`: Grid resolution for continuous distributions (default: 2048)
- `eps`: Epsilon for tail trimming (default: 1e-6)

**Returns:** Tuple `(lower, upper)`

**Examples:**
```julia
# HDI for symmetric distribution (same as equal-tailed)
d = Normal(0, 1)
hdi(d, 0.95)          # ≈ (-1.96, 1.96)

# HDI for skewed distribution (narrower than equal-tailed)
d = Exponential(1.0)
hdi(d, 0.95)          # Narrower interval
interval(d, 0.05)     # Wider interval (equal-tailed)

# HDI for discrete distribution
d = Binomial(10, 0.3)
hdi(d, 0.9)           # Smallest set of values containing 90%
```

---

## Comparing Interval Types

### Symmetric Distributions

For symmetric distributions (Normal, Cauchy, Logistic, etc.), both methods give the same result:

```julia
d = Normal(0, 1)

interval(d, 0.05)     # (-1.96, 1.96)
hdi(d, 0.95)          # (-1.96, 1.96) - same!
```

### Skewed Distributions

For skewed distributions, HDI is narrower:

```julia
d = Exponential(1.0)

# Equal-tailed: wider
lo_eq, hi_eq = interval(d, 0.05)
width_eq = hi_eq - lo_eq

# HDI: narrower
lo_hdi, hi_hdi = hdi(d, 0.95)
width_hdi = hi_hdi - lo_hdi

println("Equal-tailed width: $width_eq")
println("HDI width: $width_hdi")
# HDI is always ≤ equal-tailed width
```

### Visual Comparison

```julia
using UniDist

# Gamma distribution (skewed)
d = Gamma(2, 2)

# Equal-tailed 90% interval
eq_lo, eq_hi = interval(d, 0.10)
println("Equal-tailed: [$eq_lo, $eq_hi], width = $(eq_hi - eq_lo)")

# HDI 90% interval
hdi_lo, hdi_hi = hdi(d, 0.90)
println("HDI: [$hdi_lo, $hdi_hi], width = $(hdi_hi - hdi_lo)")

# The HDI is narrower and shifted toward the mode
```

---

## When to Use Each Method

### Use Equal-Tailed Intervals When:
- Distribution is symmetric
- You want equal error probability in both tails
- Communicating with audiences familiar with confidence intervals
- Computing p-values (two-tailed tests)

### Use HDI When:
- Distribution is skewed
- You want the shortest possible interval
- Doing Bayesian inference
- The mode is more relevant than the median
- Distribution may be multimodal

---

## Practical Applications

### Confidence Intervals for Parameters

```julia
# Estimating population mean
# After Bayesian analysis, posterior is Normal(5.2, 0.8)
posterior = Normal(5.2, 0.8)

# 95% credible interval
lo, hi = interval(posterior, 0.05)
println("95% CI for mean: [$lo, $hi]")
```

### Bayesian Inference

```julia
# Beta posterior for success probability
# Prior: Beta(1, 1), Data: 7 successes, 3 failures
posterior = Beta(1 + 7, 1 + 3)  # Beta(8, 4)

# 95% HDI
lo, hi = hdi(posterior, 0.95)
println("95% HDI for p: [$lo, $hi]")

# Compare with equal-tailed
lo_eq, hi_eq = interval(posterior, 0.05)
println("95% Equal-tailed: [$lo_eq, $hi_eq]")
```

### Predictive Intervals

```julia
# Predict next observation from Poisson(λ=5)
d = Poisson(5.0)

# 90% prediction interval
lo, hi = hdi(d, 0.90)
println("90% of observations will be in [$lo, $hi]")
```

### Comparing Groups

```julia
# Effect size has posterior Normal(0.5, 0.2)
effect = Normal(0.5, 0.2)

# 95% HDI
lo, hi = hdi(effect, 0.95)

# Check if interval excludes zero (evidence of effect)
if lo > 0 || hi < 0
    println("95% HDI excludes zero: evidence of effect")
else
    println("95% HDI includes zero: no clear evidence")
end
```

---

## Discrete Distribution HDI

For discrete distributions, HDI finds the smallest set of values containing the desired mass:

```julia
d = Binomial(20, 0.3)

# Find 90% HDI
lo, hi = hdi(d, 0.90)
println("90% HDI: {$lo, ..., $hi}")

# Verify the coverage
coverage = sum(pdf(d, k) for k in lo:hi)
println("Actual coverage: $(round(coverage * 100, digits=1))%")
```

**Note:** For discrete distributions, the actual coverage may exceed the requested mass because probability is concentrated at discrete points.

---

## Algorithm Details

### Continuous HDI Algorithm

1. Create a fine grid over the distribution support
2. Compute PDF at each grid point
3. Sort points by density (highest first)
4. Include points until cumulative mass ≥ target
5. Return (min, max) of included points

### Discrete HDI Algorithm

1. Enumerate support of distribution
2. Compute PMF at each point
3. Sort by probability (highest first)
4. Include points until cumulative mass ≥ target
5. Return (min, max) of included points

---

## Limitations

- **Multimodal distributions:** HDI returns a single interval (min to max of high-density points), which may include low-density regions between modes
- **Discrete distributions:** Require finite, enumerable support
- **Numerical precision:** Continuous HDI depends on grid resolution

---

## API Reference

```@docs
UniDist.interval
UniDist.hdi
```
