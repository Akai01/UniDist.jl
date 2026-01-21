# [Discrete Distributions](@id discrete-distributions)

UniDist.jl provides 19 discrete probability distributions. This page documents each distribution with its parameters, support, and usage examples.

## Common Distributions

### Bernoulli

Single trial with two outcomes (success/failure).

```julia
Bernoulli(p)  # p = probability of success
```

**Parameters:**
- `p`: Probability of success, 0 ≤ p ≤ 1

**Support:** {0, 1}

**Example:**
```julia
d = Bernoulli(0.7)

pdf(d, 1)             # 0.7 (probability of success)
pdf(d, 0)             # 0.3 (probability of failure)
cdf(d, 0)             # 0.3
mean(d)               # 0.7
var(d)                # 0.21
```

**Use cases:** Coin flips, yes/no outcomes, binary classification.

---

### Binomial

Number of successes in n independent Bernoulli trials.

```julia
Binomial(n, p)  # n = trials, p = success probability
```

**Parameters:**
- `n`: Number of trials, n ≥ 0
- `p`: Probability of success, 0 ≤ p ≤ 1

**Support:** {0, 1, 2, ..., n}

**Example:**
```julia
d = Binomial(10, 0.5)

pdf(d, 5)             # ≈ 0.246 (probability of exactly 5 successes)
cdf(d, 5)             # ≈ 0.623 (probability of at most 5 successes)
mean(d)               # 5.0
var(d)                # 2.5
```

**Use cases:** Quality control, survey sampling, A/B testing, clinical trials.

---

### Poisson

Number of events in a fixed interval when events occur at constant rate.

```julia
Poisson(λ)  # λ = rate (expected count)
```

**Parameters:**
- `λ` (lambda): Rate parameter, λ > 0

**Support:** {0, 1, 2, 3, ...}

**Example:**
```julia
d = Poisson(3.0)

pdf(d, 0)             # ≈ 0.050 (P(X = 0))
pdf(d, 3)             # ≈ 0.224 (P(X = 3))
cdf(d, 3)             # ≈ 0.647
mean(d)               # 3.0
var(d)                # 3.0
```

**Use cases:** Customer arrivals, website visits, defect counts, rare events.

---

### Geometric

Number of failures before first success.

```julia
Geometric(p)  # p = success probability
```

**Parameters:**
- `p`: Probability of success, 0 < p ≤ 1

**Support:** {0, 1, 2, 3, ...}

**Example:**
```julia
d = Geometric(0.5)

pdf(d, 0)             # 0.5 (success on first trial)
pdf(d, 1)             # 0.25 (one failure, then success)
pdf(d, 2)             # 0.125
mean(d)               # 1.0 (expected failures before success)
```

**Use cases:** Number of trials until success, waiting times, reliability testing.

---

### DiscreteUniform

All integers in [a, b] equally likely.

```julia
DiscreteUniform(a, b)  # a = min, b = max
```

**Parameters:**
- `a`: Minimum value (integer)
- `b`: Maximum value (integer), b ≥ a

**Support:** {a, a+1, ..., b}

**Example:**
```julia
d = DiscreteUniform(1, 6)  # Fair die

pdf(d, 3)             # ≈ 0.167
cdf(d, 3)             # 0.5
mean(d)               # 3.5
```

**Use cases:** Dice rolls, random selection, lottery numbers.

---

## Count Distributions

### Pascal (Negative Binomial)

Number of failures before r successes.

```julia
Pascal(r, p)  # r = successes needed, p = success probability
```

**Parameters:**
- `r`: Number of successes required, r > 0
- `p`: Probability of success, 0 < p ≤ 1

**Support:** {0, 1, 2, 3, ...}

**Example:**
```julia
d = Pascal(3, 0.5)

pdf(d, 2)             # Probability of 2 failures before 3 successes
mean(d)               # r(1-p)/p
```

**Use cases:** Overdispersed count data, number of trials until r successes.

---

### GammaPoisson (Negative Binomial alternative parameterization)

Poisson distribution with Gamma-distributed rate.

```julia
GammaPoisson(r, p)
```

**Use cases:** Modeling count data with overdispersion.

---

### Hypergeometric

Sampling without replacement from a finite population.

```julia
Hypergeometric(N, K, n)  # N = population, K = successes in population, n = draws
```

**Parameters:**
- `N`: Population size
- `K`: Number of success states in population
- `n`: Number of draws

**Support:** {max(0, n+K-N), ..., min(n, K)}

**Example:**
```julia
# Drawing 5 cards, how many aces?
d = Hypergeometric(52, 4, 5)

pdf(d, 0)             # Probability of no aces
pdf(d, 1)             # Probability of exactly 1 ace
mean(d)               # n * K / N
```

**Use cases:** Quality sampling, card games, capture-recapture studies.

---

### NegativeHypergeometric

```julia
NegativeHypergeometric(N, K, r)
```

---

## Special Distributions

### BetaBinomial

Binomial with Beta-distributed success probability.

```julia
BetaBinomial(n, α, β)  # n = trials, α, β = Beta parameters
```

**Parameters:**
- `n`: Number of trials
- `α` (alpha): Beta shape parameter, α > 0
- `β` (beta): Beta shape parameter, β > 0

**Support:** {0, 1, 2, ..., n}

**Example:**
```julia
d = BetaBinomial(10, 2, 3)

pdf(d, 5)             # Probability of 5 successes
mean(d)               # n * α / (α + β)
```

**Use cases:** Overdispersed binomial data, Bayesian inference.

---

### BetaPascal

Pascal distribution with Beta-distributed success probability.

```julia
BetaPascal(r, α, β)
```

---

### Polya (Beta-Binomial alternative)

```julia
Polya(n, α, β)
```

---

## Rank and Order Distributions

### Benford

Distribution of first digits in many real-world datasets.

```julia
Benford()
```

**Support:** {1, 2, 3, ..., 9}

**Example:**
```julia
d = Benford()

pdf(d, 1)             # ≈ 0.301 (30.1% start with 1)
pdf(d, 9)             # ≈ 0.046 (4.6% start with 9)
```

**Use cases:** Fraud detection, data validation, forensic accounting.

---

### Zipf

Power-law distribution for ranked data.

```julia
Zipf(s, N)  # s = exponent, N = number of elements
```

**Parameters:**
- `s`: Exponent parameter, s > 0
- `N`: Number of elements

**Support:** {1, 2, 3, ..., N}

**Example:**
```julia
d = Zipf(1.0, 100)

pdf(d, 1)             # Probability of rank 1
pdf(d, 2)             # Probability of rank 2 (about half of rank 1)
```

**Use cases:** Word frequencies, city populations, website popularity.

---

### Zeta (Zipf with infinite support)

```julia
Zeta(s)  # s = exponent
```

**Parameters:**
- `s`: Exponent parameter, s > 1

**Support:** {1, 2, 3, ...}

**Use cases:** Power-law phenomena with unbounded support.

---

### Logarithm (Logarithmic Series)

```julia
Logarithm(p)  # p = parameter
```

**Parameters:**
- `p`: Parameter, 0 < p < 1

**Support:** {1, 2, 3, ...}

**Use cases:** Species abundance, word frequency.

---

## Other Discrete Distributions

### DiscreteWeibull

Discrete analog of Weibull distribution.

```julia
DiscreteWeibull(q, β)  # q = parameter, β = shape
```

**Support:** {0, 1, 2, 3, ...}

**Use cases:** Discrete lifetime data, count-based reliability.

---

### PowerSeries

General power series distribution.

```julia
PowerSeries(f, θ)  # f = coefficient function, θ = parameter
```

---

### Rectangular

Alias for DiscreteUniform.

```julia
Rectangular(a, b)
```

---

## Usage Examples

### Comparing Distributions

```julia
using UniDist

# Compare Binomial approximation to Poisson
n, p = 100, 0.03
binomial = Binomial(n, p)
poisson = Poisson(n * p)

for k in 0:10
    println("k=$k: Binomial=$(round(pdf(binomial, k), digits=4)), Poisson=$(round(pdf(poisson, k), digits=4))")
end
```

### Cumulative Probabilities

```julia
d = Poisson(5.0)

# P(X ≤ 3)
cdf(d, 3)

# P(X > 3) = 1 - P(X ≤ 3)
1 - cdf(d, 3)

# P(2 ≤ X ≤ 5)
cdf(d, 5) - cdf(d, 1)
```

### Finding Quantiles

```julia
d = Binomial(20, 0.5)

# Find the 95th percentile
quantile(d, 0.95)

# Find the median
quantile(d, 0.5)
```
