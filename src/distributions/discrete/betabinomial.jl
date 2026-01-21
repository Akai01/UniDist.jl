import SpecialFunctions: loggamma

struct BetaBinomial{Ta,Tb,Tn} <: DiscreteDistribution
    a::Ta
    b::Tb
    n::Tn
end

function support(d::BetaBinomial)
    return d.n isa AbstractArray ? nothing : (0:Int(d.n))
end

function _beta_binomial_pmf(a::Real, b::Real, n::Real, x::Real)
    if !is_integer_value(x) || x < 0 || x > n
        return 0.0
    end
    k = Int(x)
    logpmf = loggamma(n + 1) - loggamma(k + 1) - loggamma(n - k + 1) +
        loggamma(k + a) + loggamma(n - k + b) + loggamma(a + b) -
        loggamma(n + a + b) - loggamma(a) - loggamma(b)
    return exp(logpmf)
end

pdf(d::BetaBinomial, x) =
    broadcast((a, b, n, x) -> _beta_binomial_pmf(a, b, n, x), d.a, d.b, d.n, x)

function cdf(d::BetaBinomial, x)
    return broadcast(
        (a, b, n, x) -> _discrete_cdf_scalar(k -> _beta_binomial_pmf(a, b, n, k), 0, Int(n), x),
        d.a,
        d.b,
        d.n,
        x,
    )
end

function quantile(d::BetaBinomial, u)
    return broadcast(
        (a, b, n, u) -> _discrete_quantile_scalar(k -> _beta_binomial_pmf(a, b, n, k), 0, Int(n), u),
        d.a,
        d.b,
        d.n,
        u,
    )
end

function rand(rng::AbstractRNG, d::BetaBinomial)
    shape = broadcast_shape(d.a, d.b, d.n)
    if isempty(shape)
        return quantile(d, rand(rng))
    end
    u = rand(rng, shape)
    return quantile(d, u)
end

function rand(rng::AbstractRNG, d::BetaBinomial, n::Integer)
    shape = broadcast_shape(d.a, d.b, d.n)
    if isempty(shape)
        u = rand(rng, n)
        return quantile(d, u)
    end
    u = rand(rng, (n, shape...))
    return quantile(d, u)
end

mean(d::BetaBinomial) = @. d.n * d.a / (d.a + d.b)

function var(d::BetaBinomial)
    return @. d.n * d.a * d.b * (d.a + d.b + d.n) / ((d.a + d.b) ^ 2 * (d.a + d.b + 1))
end
