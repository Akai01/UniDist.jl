import SpecialFunctions: loggamma

struct BetaPascal{Tn,Ta,Tb} <: DiscreteDistribution
    n::Tn
    a::Ta
    b::Tb
end

support(::BetaPascal) = nothing

function _beta_pascal_pmf(n::Real, a::Real, b::Real, x::Real)
    if !is_integer_value(x) || x < 0
        return 0.0
    end
    k = Int(x)
    logpmf = loggamma(n + k) - loggamma(k + 1) - loggamma(n) +
        loggamma(n + a) + loggamma(b + k) - loggamma(n + a + b + k) -
        (loggamma(a) + loggamma(b) - loggamma(a + b))
    return exp(logpmf)
end

pdf(d::BetaPascal, x) =
    broadcast((n, a, b, x) -> _beta_pascal_pmf(n, a, b, x), d.n, d.a, d.b, x)

function cdf(d::BetaPascal, x)
    return broadcast(
        (n, a, b, x) -> _discrete_cdf_scalar(k -> _beta_pascal_pmf(n, a, b, k), 0, nothing, x),
        d.n,
        d.a,
        d.b,
        x,
    )
end

function quantile(d::BetaPascal, u)
    return broadcast(
        (n, a, b, u) -> _discrete_quantile_scalar(k -> _beta_pascal_pmf(n, a, b, k), 0, nothing, u),
        d.n,
        d.a,
        d.b,
        u,
    )
end

function rand(rng::AbstractRNG, d::BetaPascal)
    shape = broadcast_shape(d.n, d.a, d.b)
    if isempty(shape)
        return quantile(d, rand(rng))
    end
    u = rand(rng, shape)
    return quantile(d, u)
end

function rand(rng::AbstractRNG, d::BetaPascal, n::Integer)
    shape = broadcast_shape(d.n, d.a, d.b)
    if isempty(shape)
        u = rand(rng, n)
        return quantile(d, u)
    end
    u = rand(rng, (n, shape...))
    return quantile(d, u)
end

mean(d::BetaPascal) = @. d.n * d.b / (d.a - 1)

function var(d::BetaPascal)
    return @. d.n * d.b * (d.a + d.b - 1) / ((d.a - 1) * (d.a - 2)) +
        (d.n ^ 2) * d.b * (d.a + d.b - 1) / ((d.a - 1) ^ 2 * (d.a - 2))
end
