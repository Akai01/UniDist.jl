# [Basic Usage](@id basic-usage)

This page demonstrates fundamental usage patterns for UniDist.jl.

## Creating and Using Distributions

### Example 1: Working with Normal Distribution

```julia
using UniDist

# Create a Normal distribution
d = Normal(100, 15)  # IQ scores: mean=100, sd=15

# What percentage of people have IQ > 130?
1 - cdf(d, 130)       # ≈ 0.0228 (about 2.3%)

# What percentage have IQ between 85 and 115?
cdf(d, 115) - cdf(d, 85)  # ≈ 0.683 (about 68%)

# What IQ score is at the 99th percentile?
quantile(d, 0.99)     # ≈ 134.9

# Generate 10 random IQ scores
r(d, 10)
```

### Example 2: Quality Control with Binomial

```julia
using UniDist

# A factory produces items with 2% defect rate
# In a batch of 100 items, what's the probability of finding:

d = Binomial(100, 0.02)

# Exactly 2 defects?
pdf(d, 2)             # ≈ 0.273

# At most 3 defects?
cdf(d, 3)             # ≈ 0.858

# More than 5 defects? (reject the batch)
1 - cdf(d, 5)         # ≈ 0.015

# Expected number of defects
mean(d)               # 2.0
```

### Example 3: Waiting Times with Exponential

```julia
using UniDist

# Customers arrive on average every 5 minutes
d = Exponential(5.0)

# Probability of waiting more than 10 minutes
sf(d, 10)             # ≈ 0.135

# Probability of waiting less than 2 minutes
cdf(d, 2)             # ≈ 0.330

# Median waiting time
quantile(d, 0.5)      # ≈ 3.47 minutes
```

---

## Computing Multiple Values

### Batch Probability Calculations

```julia
using UniDist

d = Normal(0, 1)

# PDF at multiple points
xs = -3:0.5:3
densities = pdf(d, collect(xs))

# Print as a table
println("x\t\tpdf(x)")
for (x, p) in zip(xs, densities)
    println("$x\t\t$(round(p, digits=4))")
end
```

### Multiple Quantiles

```julia
using UniDist

d = Normal(100, 15)

# Common percentiles
percentiles = [0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99]
values = quantile(d, percentiles)

println("Percentile\tValue")
for (p, v) in zip(percentiles, values)
    println("$(Int(p*100))%\t\t$(round(v, digits=1))")
end
```

---

## Distribution Properties

### Summary Statistics

```julia
using UniDist

distributions = [
    ("Normal(0,1)", Normal(0, 1)),
    ("Exponential(2)", Exponential(2)),
    ("Gamma(2,3)", Gamma(2, 3)),
    ("Beta(2,5)", Beta(2, 5)),
]

println("Distribution\t\tMean\tVariance\tSupport")
for (name, d) in distributions
    m = round(mean(d), digits=3)
    v = round(var(d), digits=3)
    s = support(d)
    println("$name\t\t$m\t$v\t\t$s")
end
```

---

## Probability Calculations

### Finding Probabilities

```julia
using UniDist

# Exam scores follow Normal(75, 10)
scores = Normal(75, 10)

# P(score ≥ 90) - probability of an A
prob_A = 1 - cdf(scores, 90)
println("P(A grade): $(round(prob_A * 100, digits=1))%")

# P(60 ≤ score < 70) - probability of a D
prob_D = cdf(scores, 70) - cdf(scores, 60)
println("P(D grade): $(round(prob_D * 100, digits=1))%")

# P(score < 50) - probability of failing
prob_fail = cdf(scores, 50)
println("P(Fail): $(round(prob_fail * 100, digits=1))%")
```

### Inverse Problems

```julia
using UniDist

# What score do you need to be in the top 10%?
scores = Normal(75, 10)
top_10_cutoff = quantile(scores, 0.90)
println("Top 10% cutoff: $(round(top_10_cutoff, digits=1))")

# What score is the median?
median_score = quantile(scores, 0.50)
println("Median score: $median_score")
```

---

## Working with Discrete Distributions

### Dice Rolling

```julia
using UniDist

# Fair die
die = DiscreteUniform(1, 6)

# Probability of rolling a 6
pdf(die, 6)           # ≈ 0.167

# Probability of rolling ≤ 3
cdf(die, 3)           # 0.5

# Expected value
mean(die)             # 3.5
```

### Counting Successes

```julia
using UniDist

# Flip a fair coin 10 times
flips = Binomial(10, 0.5)

# Probability of exactly 5 heads
pdf(flips, 5)         # ≈ 0.246

# Probability of at least 7 heads
1 - cdf(flips, 6)     # ≈ 0.172

# Most likely outcome
# (find k that maximizes pdf)
probs = [pdf(flips, k) for k in 0:10]
mode = argmax(probs) - 1  # 5
```

---

## Random Sampling

### Basic Sampling

```julia
using UniDist

d = Normal(0, 1)

# Single sample
x = r(d)

# Multiple samples
samples = r(d, 1000)

# Sample statistics
sample_mean = sum(samples) / length(samples)
sample_var = sum((x - sample_mean)^2 for x in samples) / (length(samples) - 1)

println("Sample mean: $(round(sample_mean, digits=3)) (true: 0)")
println("Sample var: $(round(sample_var, digits=3)) (true: 1)")
```

### Reproducible Sampling

```julia
using UniDist
using Random

d = Normal(0, 1)

# Set seed for reproducibility
rng = MersenneTwister(42)

# These will always be the same
samples1 = r(d, 5; rng=MersenneTwister(42))
samples2 = r(d, 5; rng=MersenneTwister(42))

samples1 == samples2  # true
```

---

## Combining Operations

### Monte Carlo Estimation

```julia
using UniDist

# Estimate P(X > Y) where X ~ Normal(5, 2), Y ~ Normal(4, 1)
X = Normal(5, 2)
Y = Normal(4, 1)

n = 100000
x_samples = r(X, n)
y_samples = r(Y, n)

prob_x_greater = sum(x_samples .> y_samples) / n
println("P(X > Y) ≈ $(round(prob_x_greater, digits=3))")
```

### Distribution Comparison

```julia
using UniDist

# Compare Poisson and Normal approximation
λ = 50
poisson = Poisson(λ)
normal_approx = Normal(λ, sqrt(λ))

# Compare CDFs at various points
points = [40, 45, 50, 55, 60]
println("x\tPoisson\t\tNormal")
for x in points
    p_pois = round(cdf(poisson, x), digits=4)
    p_norm = round(cdf(normal_approx, x), digits=4)
    println("$x\t$p_pois\t\t$p_norm")
end
```
