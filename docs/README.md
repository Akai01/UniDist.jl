# Distributional.jl

This package provides vectorized probability distributions with a consistent API across discrete and continuous families.

## API overview

Core functions:

- `pdf(dist, x)` / `cdf(dist, x)` / `quantile(dist, u)`
- Shorthands: `d`, `p`, `q`
- `rand(rng, dist)` / `rand(rng, dist, n)`
- `mean(dist)` / `var(dist)`
- `interval(dist, alpha=0.05)`
- `hdi(dist, mass=0.9)`

All functions are vectorized. Parameters can be scalars or arrays and use broadcast semantics.

## Numeric CDF/quantile

Some continuous distributions list CDF/quantile as mathematically intractable in the PDF references. For these, the implementation uses numerical integration (QuadGK) and bisection. A warning is emitted when numeric methods are used.

If you prefer to avoid warnings, you can wrap calls in `Logging.with_logger` or set a custom logger.

## Parameter conventions

The formulas follow the reference PDFs in the `Continuous/` and `Discrete/` folders. Parameter constraints are not enforced at construction time; invalid parameters may return `NaN`/`Inf` or nonsensical results.

## Example usage

```julia
using Distributional

# Bernoulli
p = [0.2, 0.6, 0.9]
bern = BernoulliQ(p)
mean(bern)           # vectorized
cdf(bern, [0, 1, 1])

# Normal
n = NormalQ(0.0, 1.0)
pdf(n, 0.0)
quantile(n, 0.95)

# Beta (numerical CDF/quantile)
b = BetaQ(2.0, 3.0)
cdf(b, 0.5)
```

## Tests

Run tests with:

```bash
julia --project -e 'using Pkg; Pkg.test()'
```
