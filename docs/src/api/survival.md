# [Survival Analysis](@id survival-analysis)

UniDist.jl provides built-in functions for survival analysis, commonly used in reliability engineering, medical research, and actuarial science.

## Overview

Survival analysis deals with time-to-event data. The key functions are:

| Function | Formula | Description |
|----------|---------|-------------|
| `sf(d, x)` | S(x) = 1 - F(x) | Survival function |
| `hazard(d, x)` | h(x) = f(x) / S(x) | Hazard function |
| `cumhaz(d, x)` | H(x) = -log(S(x)) | Cumulative hazard |

where F(x) is the CDF and f(x) is the PDF.

---

## Survival Function

### `sf` - Survival Function

The survival function S(x) = P(X > x) gives the probability of surviving beyond time x.

```julia
sf(dist, x)
```

**Formula:** `sf(d, x) = 1 - cdf(d, x)`

**Examples:**
```julia
# Exponential lifetime with mean 10
d = Exponential(10.0)

sf(d, 5)              # ≈ 0.607 (probability of surviving past t=5)
sf(d, 10)             # ≈ 0.368 (probability of surviving past t=10)
sf(d, 20)             # ≈ 0.135

# Verify: sf + cdf = 1
sf(d, 5) + cdf(d, 5)  # 1.0
```

**Interpretation:**
- sf(d, 0) = 1 (everyone alive at start)
- sf(d, ∞) = 0 (eventually everyone fails)
- Monotonically decreasing

---

## Hazard Function

### `hazard` - Instantaneous Failure Rate

The hazard function h(x) represents the instantaneous failure rate at time x, given survival up to time x.

```julia
hazard(dist, x)
```

**Formula:** `hazard(d, x) = pdf(d, x) / sf(d, x)`

**Examples:**
```julia
# Exponential: constant hazard (memoryless property)
d = Exponential(10.0)
hazard(d, 1)          # 0.1
hazard(d, 5)          # 0.1
hazard(d, 100)        # 0.1 (same at all times!)

# Weibull: flexible hazard shapes
# Shape < 1: decreasing hazard (infant mortality)
d1 = Weibull(0.5, 10)
hazard(d1, 1)         # High
hazard(d1, 10)        # Lower

# Shape = 1: constant hazard (same as Exponential)
d2 = Weibull(1.0, 10)
hazard(d2, 1)         # Constant
hazard(d2, 10)        # Same

# Shape > 1: increasing hazard (wear-out)
d3 = Weibull(2.0, 10)
hazard(d3, 1)         # Low
hazard(d3, 10)        # Higher
```

**Interpretation:**
- Units: failures per unit time
- Can increase, decrease, or remain constant
- Key for reliability modeling

---

## Cumulative Hazard Function

### `cumhaz` - Cumulative Hazard

The cumulative hazard H(x) is the integral of the hazard function from 0 to x.

```julia
cumhaz(dist, x)
```

**Formula:** `cumhaz(d, x) = -log(sf(d, x))`

**Examples:**
```julia
d = Exponential(10.0)

cumhaz(d, 5)          # 0.5
cumhaz(d, 10)         # 1.0
cumhaz(d, 20)         # 2.0

# Relationship with survival function
exp(-cumhaz(d, 5))    # ≈ sf(d, 5)
```

**Properties:**
- cumhaz(d, 0) = 0
- Monotonically increasing
- cumhaz(d, x) → ∞ as x → ∞

---

## Practical Applications

### Reliability Analysis

```julia
# Component lifetime follows Weibull distribution
# Shape = 2 indicates wear-out failure mode
component = Weibull(2.0, 1000)  # Scale = 1000 hours

# Probability of surviving 500 hours
sf(component, 500)    # ≈ 0.779

# Probability of failing between 500 and 1000 hours
sf(component, 500) - sf(component, 1000)  # ≈ 0.411

# Hazard rate at 500 hours
hazard(component, 500)  # Instantaneous failure rate
```

### Medical Survival Analysis

```julia
# Patient survival follows Exponential(median_survival)
# If median survival is 24 months
median_survival = 24
d = Exponential(median_survival / log(2))

# 1-year survival probability
sf(d, 12)             # Probability of surviving 1 year

# 5-year survival probability
sf(d, 60)             # Probability of surviving 5 years

# Median survival time
quantile(d, 0.5)      # Time at which 50% have survived
```

### Comparing Distributions

```julia
# Compare hazard patterns
exp_dist = Exponential(10)
weibull_increasing = Weibull(2.0, 10)
weibull_decreasing = Weibull(0.5, 10)

times = [1, 5, 10, 15, 20]

println("Time | Exponential | Weibull(2) | Weibull(0.5)")
for t in times
    h_exp = round(hazard(exp_dist, t), digits=4)
    h_inc = round(hazard(weibull_increasing, t), digits=4)
    h_dec = round(hazard(weibull_decreasing, t), digits=4)
    println("$t    | $h_exp      | $h_inc     | $h_dec")
end
```

### Calculating Expected Remaining Lifetime

```julia
# Mean residual life at time t
# E[X - t | X > t]
function mean_residual_life(d, t, upper=1000)
    # Approximate by numerical integration
    integrand(x) = sf(d, x) / sf(d, t)
    # Integration from t to upper
    # (simplified approximation)
    dx = 0.1
    xs = t:dx:upper
    sum(sf(d, x) / sf(d, t) * dx for x in xs)
end

d = Exponential(10)
mean_residual_life(d, 0)   # ≈ 10 (memoryless: same as original mean)
mean_residual_life(d, 5)   # ≈ 10 (memoryless property!)
```

---

## Vectorized Operations

All survival functions support vectorized operations:

```julia
d = Weibull(2.0, 100)
times = [10, 25, 50, 75, 100, 150, 200]

# Survival probabilities at multiple times
sf(d, times)

# Hazard rates at multiple times
hazard(d, times)

# Cumulative hazard at multiple times
cumhaz(d, times)
```

---

## Common Survival Distributions

| Distribution | Hazard Shape | Use Case |
|--------------|--------------|----------|
| Exponential | Constant | Memoryless failures |
| Weibull (β < 1) | Decreasing | Infant mortality |
| Weibull (β = 1) | Constant | Random failures |
| Weibull (β > 1) | Increasing | Wear-out failures |
| LogNormal | Non-monotonic | Fatigue failures |
| Gamma | Flexible | Repair times |
| Gompertz | Increasing | Human mortality |
