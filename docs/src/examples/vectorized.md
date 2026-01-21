# [Vectorized Operations](@id vectorized-operations)

UniDist.jl supports vectorized parameters, allowing efficient batch computations across multiple distribution configurations.

## Vectorized Distribution Parameters

### Multiple Means

```julia
using UniDist

# Create Normal distributions with different means
# All with σ = 1
d = Normal([0.0, 1.0, 2.0, 3.0], 1.0)

# Evaluate PDF at x = 1 for all distributions
pdf(d, 1.0)
# Returns: [0.242, 0.399, 0.242, 0.054]
# (density at x=1 for each mean)
```

### Multiple Standard Deviations

```julia
using UniDist

# Normal(0, σ) for different σ values
d = Normal(0.0, [0.5, 1.0, 2.0, 4.0])

# Compare densities at x = 0
pdf(d, 0.0)
# Higher density for smaller σ (more concentrated)
```

### Grid of Parameters

```julia
using UniDist

# Create a grid of Beta distributions
αs = [0.5, 1.0, 2.0, 5.0]
βs = [0.5, 1.0, 2.0, 5.0]

# Evaluate mean for each combination
for α in αs
    for β in βs
        d = Beta(α, β)
        m = mean(d)
        println("Beta($α, $β): mean = $(round(m, digits=3))")
    end
end
```

---

## Vectorized Evaluation Points

### PDF at Multiple Points

```julia
using UniDist

d = Normal(0, 1)

# Evaluate at many points simultaneously
xs = range(-4, 4, length=100)
densities = pdf(d, collect(xs))

# Find the maximum density
max_density = maximum(densities)
max_idx = argmax(densities)
println("Maximum density $(round(max_density, digits=4)) at x = $(xs[max_idx])")
```

### CDF at Multiple Points

```julia
using UniDist

d = Exponential(5.0)

# Compute survival curve
times = [0, 1, 2, 5, 10, 20, 50]
survival_probs = sf(d, times)

println("Time\tSurvival Probability")
for (t, s) in zip(times, survival_probs)
    println("$t\t$(round(s * 100, digits=1))%")
end
```

### Multiple Quantiles

```julia
using UniDist

d = Normal(100, 15)

# Compute percentile table
percentiles = 0.01:0.01:0.99
values = quantile(d, collect(percentiles))

# Find specific percentiles
p25 = quantile(d, 0.25)
p50 = quantile(d, 0.50)
p75 = quantile(d, 0.75)

println("Q1: $p25, Median: $p50, Q3: $p75")
println("IQR: $(p75 - p25)")
```

---

## Combining Vectorized Parameters and Points

### Full Grid Evaluation

```julia
using UniDist

# Different exponential rates
rates = [0.5, 1.0, 2.0]
times = [0.5, 1.0, 2.0, 5.0]

println("Survival probabilities S(t) = P(T > t)")
println("\nTime\tλ=0.5\tλ=1.0\tλ=2.0")

for t in times
    probs = [round(sf(Exponential(1/λ), t), digits=3) for λ in rates]
    println("$t\t$(probs[1])\t$(probs[2])\t$(probs[3])")
end
```

---

## Practical Applications

### Sensitivity Analysis

```julia
using UniDist

# How does changing parameters affect the 95th percentile?
base_mean = 100
base_sd = 15

# Vary mean
means = 90:5:110
for μ in means
    d = Normal(μ, base_sd)
    p95 = quantile(d, 0.95)
    println("μ=$μ: 95th percentile = $(round(p95, digits=1))")
end

println()

# Vary standard deviation
sds = 10:2:20
for σ in sds
    d = Normal(base_mean, σ)
    p95 = quantile(d, 0.95)
    println("σ=$σ: 95th percentile = $(round(p95, digits=1))")
end
```

### Comparing Distribution Families

```julia
using UniDist

# Compare different distributions with same mean and variance
# Mean = 5, Variance = 5

distributions = [
    ("Gamma", Gamma(5, 1)),           # shape=5, scale=1 → mean=5, var=5
    ("Normal", Normal(5, sqrt(5))),   # mean=5, var=5
]

x_vals = 0:0.5:15

println("x\tGamma\t\tNormal")
for x in x_vals
    vals = [round(pdf(d, x), digits=4) for (_, d) in distributions]
    println("$x\t$(vals[1])\t\t$(vals[2])")
end
```

### Batch Risk Calculations

```julia
using UniDist

# Calculate Value at Risk (VaR) for different confidence levels
# Returns follow Normal(0.05, 0.20) - 5% expected return, 20% volatility

returns = Normal(0.05, 0.20)
confidence_levels = [0.90, 0.95, 0.99]

println("VaR at different confidence levels:")
for α in confidence_levels
    var = -quantile(returns, 1 - α)
    println("$(Int(α*100))% VaR: $(round(var * 100, digits=2))%")
end
```

---

## Performance Tips

### Pre-allocating Results

```julia
using UniDist

d = Normal(0, 1)
n = 10000
xs = randn(n)

# Efficient: single vectorized call
result = pdf(d, xs)
```

### Avoiding Repeated Distribution Creation

```julia
using UniDist

# Less efficient: creating distribution in loop
function slow_version(means, x)
    [pdf(Normal(μ, 1), x) for μ in means]
end

# More efficient: vectorized parameters
function fast_version(means, x)
    d = Normal(means, 1)
    pdf(d, x)
end
```

---

## Working with Arrays

### Element-wise Operations

```julia
using UniDist

# Array of different distribution types
distributions = [Normal(0, 1), Normal(0, 2), Normal(1, 1)]

# Evaluate each at x = 0.5
results = [pdf(d, 0.5) for d in distributions]

# Or using broadcasting syntax
x = 0.5
results = pdf.(distributions, x)
```

### Matrix of Probabilities

```julia
using UniDist

# Create probability matrix
# Rows: different distributions
# Columns: different x values

dists = [Normal(0, 1), Normal(0, 2), Exponential(1)]
xs = [0.0, 0.5, 1.0, 2.0]

prob_matrix = [pdf(d, x) for d in dists, x in xs]

# Display as table
println("Distribution\t", join(xs, "\t"))
for (i, d) in enumerate(["N(0,1)", "N(0,2)", "Exp(1)"])
    row = [round(prob_matrix[i, j], digits=3) for j in 1:length(xs)]
    println("$d\t\t", join(row, "\t"))
end
```
