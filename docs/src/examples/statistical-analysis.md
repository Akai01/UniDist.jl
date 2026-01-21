# [Statistical Analysis](@id statistical-analysis)

This page demonstrates using UniDist.jl for common statistical analysis tasks.

## Hypothesis Testing

### One-Sample Z-Test

```julia
using UniDist

# Test if sample mean differs from hypothesized value
# H₀: μ = 100, H₁: μ ≠ 100

sample_mean = 103.5
sample_size = 50
known_sd = 15
hypothesized_mean = 100

# Calculate z-statistic
se = known_sd / sqrt(sample_size)
z = (sample_mean - hypothesized_mean) / se

# P-value (two-tailed)
standard_normal = Normal(0, 1)
p_value = 2 * (1 - cdf(standard_normal, abs(z)))

println("Z-statistic: $(round(z, digits=3))")
println("P-value: $(round(p_value, digits=4))")

if p_value < 0.05
    println("Reject H₀ at α = 0.05")
else
    println("Fail to reject H₀ at α = 0.05")
end
```

### Chi-Square Goodness of Fit

```julia
using UniDist

# Test if observed frequencies match expected
observed = [45, 35, 20]  # Observed counts
expected = [40, 40, 20]  # Expected counts

# Chi-square statistic
chi_sq = sum((o - e)^2 / e for (o, e) in zip(observed, expected))

# Degrees of freedom = categories - 1
df = length(observed) - 1

# P-value
chi_dist = ChiSquare(df)
p_value = 1 - cdf(chi_dist, chi_sq)

println("χ² statistic: $(round(chi_sq, digits=3))")
println("Degrees of freedom: $df")
println("P-value: $(round(p_value, digits=4))")
```

### F-Test for Variance Ratio

```julia
using UniDist

# Test if two populations have equal variances
var1, n1 = 25.0, 30  # Sample 1: variance, size
var2, n2 = 16.0, 25  # Sample 2: variance, size

# F-statistic (larger variance in numerator)
f_stat = var1 / var2

# Degrees of freedom
df1 = n1 - 1
df2 = n2 - 1

# P-value (two-tailed)
f_dist = F(df1, df2)
p_value = 2 * min(cdf(f_dist, f_stat), 1 - cdf(f_dist, f_stat))

println("F-statistic: $(round(f_stat, digits=3))")
println("P-value: $(round(p_value, digits=4))")
```

---

## Confidence Intervals

### CI for Population Mean (Known Variance)

```julia
using UniDist

sample_mean = 75.3
known_sd = 10.0
n = 40
confidence = 0.95

# Z critical value
z_crit = quantile(Normal(0, 1), 1 - (1 - confidence) / 2)

# Margin of error
se = known_sd / sqrt(n)
margin = z_crit * se

# Confidence interval
ci_lower = sample_mean - margin
ci_upper = sample_mean + margin

println("$(Int(confidence*100))% CI: [$(round(ci_lower, digits=2)), $(round(ci_upper, digits=2))]")
```

### CI for Population Proportion

```julia
using UniDist

# Survey: 120 out of 400 prefer option A
successes = 120
n = 400
confidence = 0.95

# Sample proportion
p_hat = successes / n

# Standard error
se = sqrt(p_hat * (1 - p_hat) / n)

# Z critical value
z_crit = quantile(Normal(0, 1), 1 - (1 - confidence) / 2)

# Confidence interval
ci_lower = p_hat - z_crit * se
ci_upper = p_hat + z_crit * se

println("Sample proportion: $(round(p_hat, digits=3))")
println("$(Int(confidence*100))% CI: [$(round(ci_lower, digits=3)), $(round(ci_upper, digits=3))]")
```

---

## Power Analysis

### Sample Size for Desired Power

```julia
using UniDist

# Detect effect size d = 0.5 with 80% power at α = 0.05

effect_size = 0.5
alpha = 0.05
power = 0.80

# Critical values
z_alpha = quantile(Normal(0, 1), 1 - alpha/2)  # Two-tailed
z_beta = quantile(Normal(0, 1), power)

# Required sample size per group (two-sample t-test approximation)
n = 2 * ((z_alpha + z_beta) / effect_size)^2

println("Required sample size per group: $(ceil(Int, n))")
```

### Power for Given Sample Size

```julia
using UniDist

# Calculate power for n=50 per group, effect size=0.4, α=0.05

n = 50
effect_size = 0.4
alpha = 0.05

# Non-centrality parameter
ncp = effect_size * sqrt(n / 2)

# Critical value under null
z_crit = quantile(Normal(0, 1), 1 - alpha/2)

# Power = P(reject H₀ | H₁ true)
# Under alternative, test statistic ~ N(ncp, 1)
alt_dist = Normal(ncp, 1)
power = 1 - cdf(alt_dist, z_crit) + cdf(alt_dist, -z_crit)

println("Power: $(round(power * 100, digits=1))%")
```

---

## Distribution Fitting Assessment

### Q-Q Plot Data

```julia
using UniDist

# Generate theoretical quantiles for Q-Q plot
data = [2.3, 3.1, 3.5, 4.2, 4.8, 5.1, 5.9, 6.3, 7.1, 8.2]
n = length(data)

# Sort data
sorted_data = sort(data)

# Theoretical quantiles (standard normal)
theoretical_dist = Normal(0, 1)
probabilities = [(i - 0.5) / n for i in 1:n]
theoretical_quantiles = quantile(theoretical_dist, probabilities)

println("Data Quantile\tTheoretical Quantile")
for (d, t) in zip(sorted_data, theoretical_quantiles)
    println("$(round(d, digits=2))\t\t$(round(t, digits=2))")
end
```

### Probability Plot Correlation

```julia
using UniDist

# Check if data follows exponential distribution
data = [0.5, 0.8, 1.2, 1.5, 2.1, 2.8, 3.5, 4.2, 5.1, 7.3]
n = length(data)

# Fit exponential: estimate rate from data
mean_data = sum(data) / n
fitted_dist = Exponential(mean_data)

# Theoretical quantiles
sorted_data = sort(data)
probs = [(i - 0.5) / n for i in 1:n]
theoretical = quantile(fitted_dist, probs)

# Correlation coefficient
mean_sorted = sum(sorted_data) / n
mean_theoretical = sum(theoretical) / n

numerator = sum((sorted_data[i] - mean_sorted) * (theoretical[i] - mean_theoretical) for i in 1:n)
denom_sorted = sqrt(sum((x - mean_sorted)^2 for x in sorted_data))
denom_theoretical = sqrt(sum((x - mean_theoretical)^2 for x in theoretical))

correlation = numerator / (denom_sorted * denom_theoretical)

println("Probability plot correlation: $(round(correlation, digits=4))")
println("(Values close to 1 suggest good fit)")
```

---

## Risk Analysis

### Value at Risk (VaR)

```julia
using UniDist

# Portfolio returns assumed to follow Normal distribution
# Mean daily return: 0.05%, Volatility: 2%
daily_return = Normal(0.0005, 0.02)

# VaR at different confidence levels
println("Value at Risk (daily):")
for conf in [0.90, 0.95, 0.99]
    var = -quantile(daily_return, 1 - conf)
    println("$(Int(conf*100))% VaR: $(round(var * 100, digits=2))%")
end
```

### Expected Shortfall (CVaR)

```julia
using UniDist

# For Normal distribution, ES has closed form
μ = 0.0005
σ = 0.02
returns = Normal(μ, σ)

alpha = 0.05  # 95% confidence

# VaR
var = -quantile(returns, alpha)

# Expected Shortfall for Normal: ES = μ + σ * φ(Φ⁻¹(α)) / α
# where φ is PDF and Φ⁻¹ is quantile of standard normal
z_alpha = quantile(Normal(0, 1), alpha)
es = -(μ - σ * pdf(Normal(0, 1), z_alpha) / alpha)

println("95% VaR: $(round(var * 100, digits=2))%")
println("95% ES (CVaR): $(round(es * 100, digits=2))%")
```

---

## Quality Control

### Process Capability

```julia
using UniDist

# Process specifications
lsl = 95   # Lower specification limit
usl = 105  # Upper specification limit

# Process parameters (from sample)
process_mean = 100.2
process_sd = 1.5

process = Normal(process_mean, process_sd)

# Defect rates
p_below_lsl = cdf(process, lsl)
p_above_usl = 1 - cdf(process, usl)
total_defect_rate = p_below_lsl + p_above_usl

println("P(below LSL): $(round(p_below_lsl * 1e6, digits=1)) ppm")
println("P(above USL): $(round(p_above_usl * 1e6, digits=1)) ppm")
println("Total defect rate: $(round(total_defect_rate * 1e6, digits=1)) ppm")

# Process capability indices
cp = (usl - lsl) / (6 * process_sd)
cpk = min(usl - process_mean, process_mean - lsl) / (3 * process_sd)

println("\nCp: $(round(cp, digits=3))")
println("Cpk: $(round(cpk, digits=3))")
```

### Control Chart Limits

```julia
using UniDist

# X-bar chart for sample means
# Process: μ = 50, σ = 5, sample size n = 4

mu = 50
sigma = 5
n = 4

# Distribution of sample means
se = sigma / sqrt(n)
xbar_dist = Normal(mu, se)

# 3-sigma control limits
lcl = quantile(xbar_dist, 0.00135)  # ≈ μ - 3σ/√n
ucl = quantile(xbar_dist, 0.99865)  # ≈ μ + 3σ/√n

println("Center Line: $mu")
println("LCL: $(round(lcl, digits=2))")
println("UCL: $(round(ucl, digits=2))")

# Probability of false alarm (Type I error)
p_false_alarm = 2 * cdf(xbar_dist, lcl)
println("P(false alarm): $(round(p_false_alarm * 100, digits=3))%")
```

---

## Acceptance Sampling

### Single Sampling Plan

```julia
using UniDist

# Lot size N = 1000, sample size n = 50, accept if defects ≤ c = 2
# What is the probability of accepting a lot with 5% defective?

n = 50
c = 2
p_defective = 0.05

# Number of defectives in sample follows Binomial(n, p)
defectives = Binomial(n, p_defective)

# Probability of acceptance
p_accept = cdf(defectives, c)

println("Sampling plan: n=$n, c=$c")
println("If lot is 5% defective:")
println("P(accept) = $(round(p_accept * 100, digits=1))%")
println("P(reject) = $(round((1 - p_accept) * 100, digits=1))%")
```

### Operating Characteristic Curve

```julia
using UniDist

# OC curve for sampling plan n=50, c=2
n = 50
c = 2

println("Lot defect rate\tP(accept)")
for p in [0.01, 0.02, 0.03, 0.05, 0.07, 0.10, 0.15]
    defectives = Binomial(n, p)
    p_accept = cdf(defectives, c)
    println("$(Int(p*100))%\t\t$(round(p_accept * 100, digits=1))%")
end
```
