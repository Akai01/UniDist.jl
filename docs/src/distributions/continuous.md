# [Continuous Distributions](@id continuous-distributions)

UniDist.jl provides 54 continuous probability distributions. This page documents each distribution with its parameters, support, and usage examples.

## Common Distributions

### Normal (Gaussian)

The most widely used continuous distribution, characterized by its bell-shaped curve.

```julia
Normal(μ, σ)  # μ = mean, σ = standard deviation
```

**Parameters:**
- `μ` (mu): Mean (location parameter)
- `σ` (sigma): Standard deviation (scale parameter), σ > 0

**Support:** (-∞, +∞)

**Example:**
```julia
d = Normal(0.0, 1.0)  # Standard normal

pdf(d, 0.0)           # ≈ 0.3989
cdf(d, 1.96)          # ≈ 0.975
quantile(d, 0.975)    # ≈ 1.96
mean(d)               # 0.0
var(d)                # 1.0
```

**Use cases:** Modeling measurement errors, natural phenomena, test scores, financial returns.

---

### Exponential

Models the time between events in a Poisson process.

```julia
Exponential(θ)  # θ = scale parameter (mean)
```

**Parameters:**
- `θ` (theta): Scale parameter (mean), θ > 0

**Support:** [0, +∞)

**Example:**
```julia
d = Exponential(2.0)  # Mean time = 2

pdf(d, 1.0)           # ≈ 0.303
cdf(d, 2.0)           # ≈ 0.632
mean(d)               # 2.0
```

**Use cases:** Waiting times, survival analysis, reliability engineering, radioactive decay.

---

### Uniform

All values in an interval are equally likely.

```julia
Uniform(a, b)  # a = lower bound, b = upper bound
```

**Parameters:**
- `a`: Lower bound
- `b`: Upper bound, b > a

**Support:** [a, b]

**Example:**
```julia
d = Uniform(0.0, 10.0)

pdf(d, 5.0)           # 0.1
cdf(d, 5.0)           # 0.5
quantile(d, 0.5)      # 5.0
mean(d)               # 5.0
```

**Use cases:** Random number generation, modeling complete uncertainty within bounds.

---

### Beta

Flexible distribution on [0, 1], often used for probabilities.

```julia
Beta(α, β)  # α, β = shape parameters
```

**Parameters:**
- `α` (alpha): Shape parameter, α > 0
- `β` (beta): Shape parameter, β > 0

**Support:** [0, 1]

**Example:**
```julia
d = Beta(2.0, 5.0)

pdf(d, 0.3)           # ≈ 2.06
cdf(d, 0.5)           # ≈ 0.89
mean(d)               # ≈ 0.286
```

**Use cases:** Bayesian inference for proportions, modeling probabilities, A/B testing.

---

### Gamma

Generalization of exponential distribution, models waiting times for multiple events.

```julia
Gamma(α, θ)  # α = shape, θ = scale
```

**Parameters:**
- `α` (alpha): Shape parameter, α > 0
- `θ` (theta): Scale parameter, θ > 0

**Support:** [0, +∞)

**Example:**
```julia
d = Gamma(2.0, 2.0)

pdf(d, 2.0)           # ≈ 0.184
cdf(d, 4.0)           # ≈ 0.594
mean(d)               # 4.0 (α * θ)
```

**Use cases:** Waiting times, insurance claims, rainfall amounts, Bayesian inference.

---

## Location-Scale Distributions

### Cauchy

Heavy-tailed distribution with no defined mean or variance.

```julia
Cauchy(x₀, γ)  # x₀ = location, γ = scale
```

**Parameters:**
- `x₀`: Location parameter
- `γ` (gamma): Scale parameter, γ > 0

**Support:** (-∞, +∞)

**Example:**
```julia
d = Cauchy(0.0, 1.0)  # Standard Cauchy

pdf(d, 0.0)           # ≈ 0.318
cdf(d, 0.0)           # 0.5
```

**Use cases:** Modeling outliers, ratio of normal variables, resonance in physics.

---

### Laplace (Double Exponential)

Symmetric distribution with heavier tails than normal.

```julia
Laplace(μ, b)  # μ = location, b = scale
```

**Parameters:**
- `μ` (mu): Location parameter
- `b`: Scale parameter, b > 0

**Support:** (-∞, +∞)

**Example:**
```julia
d = Laplace(0.0, 1.0)

pdf(d, 0.0)           # 0.5
cdf(d, 0.0)           # 0.5
```

**Use cases:** Finance, signal processing, LASSO regression, robust estimation.

---

### Logistic

Similar to normal but with heavier tails.

```julia
Logistic(μ, s)  # μ = location, s = scale
```

**Parameters:**
- `μ` (mu): Location (mean)
- `s`: Scale parameter, s > 0

**Support:** (-∞, +∞)

**Example:**
```julia
d = Logistic(0.0, 1.0)

pdf(d, 0.0)           # 0.25
cdf(d, 0.0)           # 0.5
```

**Use cases:** Logistic regression, growth models, neural networks.

---

## Chi-Square Family

### ChiSquare

Sum of squared standard normal variables.

```julia
ChiSquare(ν)  # ν = degrees of freedom
```

**Parameters:**
- `ν` (nu): Degrees of freedom, ν > 0

**Support:** [0, +∞)

**Example:**
```julia
d = ChiSquare(5)

pdf(d, 3.0)           # ≈ 0.154
cdf(d, 5.0)           # ≈ 0.584
mean(d)               # 5.0
```

**Use cases:** Hypothesis testing, confidence intervals, goodness-of-fit tests.

---

### Chi

Square root of chi-square distribution.

```julia
Chi(ν)  # ν = degrees of freedom
```

**Parameters:**
- `ν` (nu): Degrees of freedom, ν > 0

**Support:** [0, +∞)

---

### F

Ratio of two chi-square distributions.

```julia
F(ν₁, ν₂)  # ν₁, ν₂ = degrees of freedom
```

**Parameters:**
- `ν₁`: Numerator degrees of freedom
- `ν₂`: Denominator degrees of freedom

**Support:** [0, +∞)

**Example:**
```julia
d = F(5, 10)

pdf(d, 1.0)           # ≈ 0.61
cdf(d, 2.0)           # ≈ 0.84
```

**Use cases:** ANOVA, comparing variances, regression analysis.

---

## Extreme Value Distributions

### ExtremeValue (Gumbel)

Models the maximum of many samples.

```julia
ExtremeValue(μ, σ)  # μ = location, σ = scale
```

**Parameters:**
- `μ` (mu): Location parameter
- `σ` (sigma): Scale parameter, σ > 0

**Support:** (-∞, +∞)

**Use cases:** Extreme weather events, flood analysis, material strength.

---

### Weibull

Flexible distribution for reliability and survival analysis.

```julia
Weibull(α, θ)  # α = shape, θ = scale
```

**Parameters:**
- `α` (alpha): Shape parameter, α > 0
- `θ` (theta): Scale parameter, θ > 0

**Support:** [0, +∞)

**Example:**
```julia
d = Weibull(2.0, 1.0)

pdf(d, 0.5)           # ≈ 0.779
cdf(d, 1.0)           # ≈ 0.632
```

**Use cases:** Reliability engineering, failure analysis, wind speed modeling.

---

### Pareto

Heavy-tailed distribution following power law.

```julia
Pareto(α, θ)  # α = shape, θ = scale (minimum)
```

**Parameters:**
- `α` (alpha): Shape parameter, α > 0
- `θ` (theta): Scale (minimum value), θ > 0

**Support:** [θ, +∞)

**Use cases:** Income distribution, city populations, file sizes, insurance claims.

---

## Log-Transformed Distributions

### LogNormal

Distribution of exp(X) where X is normal.

```julia
LogNormal(μ, σ)  # μ, σ = parameters of log(X)
```

**Parameters:**
- `μ` (mu): Mean of log(X)
- `σ` (sigma): Standard deviation of log(X), σ > 0

**Support:** (0, +∞)

**Example:**
```julia
d = LogNormal(0.0, 1.0)

pdf(d, 1.0)           # ≈ 0.399
cdf(d, 1.0)           # 0.5
```

**Use cases:** Stock prices, income, biological measurements, particle sizes.

---

### LogLogistic

Log-transformed logistic distribution.

```julia
LogLogistic(α, β)  # α = scale, β = shape
```

**Support:** (0, +∞)

**Use cases:** Survival analysis, hydrology, economics.

---

## Other Continuous Distributions

### Arcsin

U-shaped distribution on [0, 1].

```julia
Arcsin(a, b)  # a = lower bound, b = upper bound
```

**Support:** [a, b]

---

### Arctangent

```julia
Arctangent(θ, σ)
```

---

### Erlang

Special case of Gamma with integer shape.

```julia
Erlang(k, θ)  # k = shape (integer), θ = scale
```

**Use cases:** Queuing theory, telecommunications.

---

### Error

```julia
Error(μ, σ, p)
```

---

### ExponentialPower (Generalized Normal)

```julia
ExponentialPower(μ, σ, p)
```

---

### GeneralizedGamma

```julia
GeneralizedGamma(a, d, p)
```

---

### GeneralizedPareto

```julia
GeneralizedPareto(μ, σ, ξ)
```

**Use cases:** Extreme value analysis, tail risk modeling.

---

### Gompertz

```julia
Gompertz(η, b)
```

**Use cases:** Mortality modeling, actuarial science.

---

### HyperbolicSecant

```julia
HyperbolicSecant(μ, σ)
```

---

### Hyperexponential

```julia
Hyperexponential(probs, rates)
```

---

### Hypoexponential

```julia
Hypoexponential(rates)
```

---

### IDB (Inverse Beta Distribution)

```julia
IDB(α, β, θ)
```

---

### InverseGaussian

```julia
InverseGaussian(μ, λ)
```

**Use cases:** First passage times, reliability.

---

### InvertedBeta

```julia
InvertedBeta(α, β)
```

---

### InvertedGamma

```julia
InvertedGamma(α, θ)
```

**Use cases:** Bayesian inference for variance.

---

### KolmogorovSmirnov

```julia
KolmogorovSmirnov()
```

**Use cases:** Goodness-of-fit testing.

---

### LogGamma

```julia
LogGamma(α, β)
```

---

### LogisticExponential

```julia
LogisticExponential(λ, κ)
```

---

### Lomax (Pareto Type II)

```julia
Lomax(α, λ)
```

---

### Makeham

```julia
Makeham(a, b, c)
```

**Use cases:** Mortality modeling.

---

### Minimax

```julia
Minimax(β, γ)
```

---

### Muth

```julia
Muth(α)
```

---

### Power

```julia
Power(α, a, b)
```

---

### Rayleigh

```julia
Rayleigh(σ)
```

**Use cases:** Wind speed, wave height, signal processing.

---

### Triangular

```julia
Triangular(a, b, c)  # a = min, b = max, c = mode
```

**Use cases:** Project management (PERT), when only min/max/mode are known.

---

### VonMises

Circular distribution for angular data.

```julia
VonMises(μ, κ)  # μ = mean direction, κ = concentration
```

**Use cases:** Wind directions, compass bearings, time-of-day data.

---

## Noncentral Distributions

### NoncentralBeta

```julia
NoncentralBeta(α, β, λ)
```

---

### NoncentralChiSquare

```julia
NoncentralChiSquare(ν, λ)
```

**Use cases:** Power analysis, signal detection.

---

### NoncentralF

```julia
NoncentralF(ν₁, ν₂, λ)
```

---

### NoncentralT

```julia
NoncentralT(ν, λ)
```

**Use cases:** Power analysis, effect size calculations.

---

### DoublyNoncentralF

```julia
DoublyNoncentralF(ν₁, ν₂, λ₁, λ₂)
```

---

### DoublyNoncentralT

```julia
DoublyNoncentralT(ν, λ₁, λ₂)
```

---

## Standard (Parameterless) Distributions

These distributions have fixed parameters for convenience:

```julia
StandardNormal()      # Normal(0, 1)
StandardUniform()     # Uniform(0, 1)
StandardCauchy()      # Cauchy(0, 1)
StandardTriangular()  # Triangular(-1, 1, 0)
StandardPower(α)      # Power distribution on [0, 1]
```

**Example:**
```julia
d = StandardNormal()
pdf(d, 0.0)           # ≈ 0.3989
cdf(d, 1.96)          # ≈ 0.975
```
