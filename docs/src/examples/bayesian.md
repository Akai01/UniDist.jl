# [Bayesian Inference](@id bayesian-inference)

This page demonstrates using UniDist.jl for Bayesian statistical analysis with conjugate priors.

## Introduction to Bayesian Analysis

Bayesian inference combines prior beliefs with observed data to produce posterior distributions:

```
Posterior ∝ Likelihood × Prior
```

UniDist.jl provides the distributions needed for common conjugate analysis.

---

## Beta-Binomial Model

The Beta distribution is conjugate to the Binomial likelihood, making it ideal for estimating proportions.

### Estimating Success Probability

```julia
using UniDist

# Scenario: Estimating website conversion rate
# Prior belief: ~10% conversion, somewhat uncertain

# Prior: Beta(2, 18) → mean = 2/20 = 0.10
prior = Beta(2, 18)

println("Prior:")
println("  Mean: $(round(mean(prior), digits=3))")
lo, hi = hdi(prior, 0.95)
println("  95% HDI: [$(round(lo, digits=3)), $(round(hi, digits=3))]")

# Data: 15 conversions out of 100 visitors
successes = 15
failures = 85

# Posterior: Beta(α + successes, β + failures)
posterior = Beta(2 + successes, 18 + failures)

println("\nPosterior (after 100 visitors, 15 conversions):")
println("  Mean: $(round(mean(posterior), digits=3))")
lo, hi = hdi(posterior, 0.95)
println("  95% HDI: [$(round(lo, digits=3)), $(round(hi, digits=3))]")

# Probability conversion rate > 10%
p_above_10 = 1 - cdf(posterior, 0.10)
println("  P(rate > 10%): $(round(p_above_10 * 100, digits=1))%")
```

### A/B Testing

```julia
using UniDist

# Compare two versions: A (control) vs B (treatment)
# Prior: Beta(1, 1) = Uniform (non-informative)

# Version A: 45 conversions out of 500
post_A = Beta(1 + 45, 1 + 455)

# Version B: 62 conversions out of 500
post_B = Beta(1 + 62, 1 + 438)

println("Version A conversion rate:")
println("  Mean: $(round(mean(post_A) * 100, digits=2))%")
lo, hi = hdi(post_A, 0.95)
println("  95% HDI: [$(round(lo * 100, digits=2))%, $(round(hi * 100, digits=2))%]")

println("\nVersion B conversion rate:")
println("  Mean: $(round(mean(post_B) * 100, digits=2))%")
lo, hi = hdi(post_B, 0.95)
println("  95% HDI: [$(round(lo * 100, digits=2))%, $(round(hi * 100, digits=2))%]")

# Monte Carlo: P(B > A)
n_samples = 100000
samples_A = r(post_A, n_samples)
samples_B = r(post_B, n_samples)
p_B_better = sum(samples_B .> samples_A) / n_samples

println("\nP(B > A): $(round(p_B_better * 100, digits=1))%")

# Expected lift
lift_samples = (samples_B .- samples_A) ./ samples_A
mean_lift = sum(lift_samples) / n_samples
println("Expected lift: $(round(mean_lift * 100, digits=1))%")
```

---

## Gamma-Poisson Model

The Gamma distribution is conjugate to the Poisson likelihood, useful for estimating rates.

### Estimating Event Rate

```julia
using UniDist

# Scenario: Estimating customer support tickets per day
# Prior: Gamma(2, 10) → mean = 20 tickets/day, moderate uncertainty

prior = Gamma(2, 10)

println("Prior:")
println("  Mean: $(round(mean(prior), digits=1)) tickets/day")
lo, hi = hdi(prior, 0.95)
println("  95% HDI: [$(round(lo, digits=1)), $(round(hi, digits=1))]")

# Data: Observed 156 tickets over 7 days
total_count = 156
n_days = 7

# Posterior: Gamma(α + total_count, β / (1 + n_days * β))
# Or equivalently with rate parameterization
α_post = 2 + total_count
θ_post = 10 / (1 + n_days * 0.1)  # Adjust scale

# Simpler: Gamma(α + Σx, θ/(1 + n*θ)) where θ is scale
posterior = Gamma(α_post, θ_post)

println("\nPosterior (after observing $total_count tickets in $n_days days):")
println("  Mean: $(round(mean(posterior), digits=1)) tickets/day")
lo, hi = hdi(posterior, 0.95)
println("  95% HDI: [$(round(lo, digits=1)), $(round(hi, digits=1))]")
```

---

## Normal-Normal Model

With known variance, the Normal distribution is conjugate to itself.

### Estimating Population Mean

```julia
using UniDist

# Scenario: Estimating average customer spending
# Prior belief: $50 average, uncertain (σ_prior = $10)
# Known population σ = $20

prior_mean = 50.0
prior_sd = 10.0
known_sd = 20.0

prior = Normal(prior_mean, prior_sd)

println("Prior:")
println("  Mean: \$$(prior_mean)")
lo, hi = interval(prior, 0.05)
println("  95% CI: [\$$(round(lo, digits=2)), \$$(round(hi, digits=2))]")

# Data: n=25 customers, sample mean = $58
n = 25
sample_mean = 58.0

# Posterior parameters
# Precision-weighted combination
prior_precision = 1 / prior_sd^2
data_precision = n / known_sd^2
total_precision = prior_precision + data_precision

posterior_mean = (prior_precision * prior_mean + data_precision * sample_mean) / total_precision
posterior_sd = sqrt(1 / total_precision)

posterior = Normal(posterior_mean, posterior_sd)

println("\nPosterior (after n=$n observations):")
println("  Mean: \$$(round(posterior_mean, digits=2))")
lo, hi = interval(posterior, 0.05)
println("  95% CI: [\$$(round(lo, digits=2)), \$$(round(hi, digits=2))]")

# Probability mean > $55
p_above_55 = 1 - cdf(posterior, 55)
println("  P(μ > \$55): $(round(p_above_55 * 100, digits=1))%")
```

---

## Sequential Updating

Bayesian inference naturally supports sequential data processing.

```julia
using UniDist

# Start with weakly informative prior
prior = Beta(1, 1)

println("Sequential Bayesian updating for proportion:")
println("=" ^ 50)

# Observe data in batches
batches = [(5, 10), (8, 15), (12, 20), (7, 12)]  # (successes, trials)

current = prior
total_successes = 0
total_trials = 0

println("\nBatch\tData\t\tPosterior Mean\t95% HDI")
println("-" ^ 50)

for (i, (s, n)) in enumerate(batches)
    total_successes += s
    total_trials += n

    # Update posterior
    current = Beta(1 + total_successes, 1 + total_trials - total_successes)

    lo, hi = hdi(current, 0.95)
    println("$i\t$s/$n\t\t$(round(mean(current), digits=3))\t\t[$(round(lo, digits=3)), $(round(hi, digits=3))]")
end

println("-" ^ 50)
println("Final: $(total_successes)/$(total_trials) total")
```

---

## Predictive Distributions

### Posterior Predictive for Binomial

```julia
using UniDist

# After observing data, predict future observations
# Posterior: Beta(20, 80) for success probability

posterior = Beta(20, 80)

# Predict: In next 50 trials, how many successes?
# This is Beta-Binomial

n_future = 50

# Monte Carlo prediction
n_samples = 10000
predictions = Int[]

for _ in 1:n_samples
    # Sample p from posterior
    p = r(posterior)[1]
    # Sample outcome from Binomial(n_future, p)
    k = r(Binomial(n_future, p))[1]
    push!(predictions, k)
end

mean_pred = sum(predictions) / n_samples
lo = sort(predictions)[Int(0.025 * n_samples)]
hi = sort(predictions)[Int(0.975 * n_samples)]

println("Posterior predictive for next $n_future trials:")
println("  Expected successes: $(round(mean_pred, digits=1))")
println("  95% prediction interval: [$lo, $hi]")
```

---

## Model Comparison

### Bayes Factor Approximation

```julia
using UniDist

# Compare two models for coin fairness
# M1: Fair coin (p = 0.5)
# M2: Biased coin (p ~ Beta(1, 1))

# Data: 7 heads out of 10 flips

k = 7  # successes
n = 10 # trials

# Likelihood under M1 (point hypothesis p = 0.5)
L1 = pdf(Binomial(n, 0.5), k)

# Marginal likelihood under M2 (integrated over prior)
# For Beta(1,1) prior and Binomial likelihood:
# P(data | M2) = B(k+1, n-k+1) / B(1, 1) = 1/(n+1)
L2 = 1 / (n + 1)

# Bayes Factor
BF_12 = L1 / L2

println("Data: $k heads out of $n flips")
println("\nP(data | fair coin): $(round(L1, digits=4))")
println("P(data | biased coin): $(round(L2, digits=4))")
println("\nBayes Factor (fair vs biased): $(round(BF_12, digits=2))")

if BF_12 > 3
    println("Evidence favors fair coin")
elseif BF_12 < 1/3
    println("Evidence favors biased coin")
else
    println("Evidence is inconclusive")
end
```

---

## Credible Intervals vs Confidence Intervals

```julia
using UniDist

# Bayesian credible interval
# Posterior for proportion after observing 30/100

posterior = Beta(31, 71)  # Beta(1+30, 1+70) with uniform prior

# 95% Equal-tailed credible interval
eq_lo, eq_hi = interval(posterior, 0.05)

# 95% Highest Density Interval
hdi_lo, hdi_hi = hdi(posterior, 0.95)

println("Posterior: Beta(31, 71)")
println("  Mean: $(round(mean(posterior), digits=3))")
println()
println("95% Equal-tailed CI: [$(round(eq_lo, digits=3)), $(round(eq_hi, digits=3))]")
println("95% HDI: [$(round(hdi_lo, digits=3)), $(round(hdi_hi, digits=3))]")
println()
println("Interpretation:")
println("  There is a 95% probability that the true proportion")
println("  lies within the credible interval, given the data and prior.")
```

---

## ROPE Analysis (Region of Practical Equivalence)

```julia
using UniDist

# Is the effect practically equivalent to zero?
# Define ROPE as [-0.1, 0.1]

posterior = Normal(0.08, 0.05)  # Effect size posterior

rope_lo = -0.1
rope_hi = 0.1

# Probability mass in ROPE
p_in_rope = cdf(posterior, rope_hi) - cdf(posterior, rope_lo)

# Probability below ROPE (negative effect)
p_below = cdf(posterior, rope_lo)

# Probability above ROPE (positive effect)
p_above = 1 - cdf(posterior, rope_hi)

println("Posterior: Normal(0.08, 0.05)")
println("ROPE: [$rope_lo, $rope_hi]")
println()
println("P(effect in ROPE): $(round(p_in_rope * 100, digits=1))%")
println("P(effect < ROPE): $(round(p_below * 100, digits=1))%")
println("P(effect > ROPE): $(round(p_above * 100, digits=1))%")
println()

if p_in_rope > 0.95
    println("Decision: Accept practical equivalence")
elseif p_above > 0.95 || p_below > 0.95
    println("Decision: Reject practical equivalence")
else
    println("Decision: Undecided (need more data)")
end
```
