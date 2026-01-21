import SpecialFunctions: loggamma

struct Pascal{Tn,Tp} <: DiscreteDistribution
    n::Tn
    p::Tp
end

support(::Pascal) = nothing

function _pascal_pmf(n::Real, p::Real, x::Real)
    if !is_integer_value(x) || x < 0
        return 0.0
    end
    k = Int(x)
    logpmf = loggamma(n + k) - loggamma(k + 1) - loggamma(n) + n * log(p) + k * log(1 - p)
    return exp(logpmf)
end

pdf(d::Pascal, x) = broadcast((n, p, x) -> _pascal_pmf(n, p, x), d.n, d.p, x)

function cdf(d::Pascal, x)
    return broadcast(
        (n, p, x) -> _discrete_cdf_scalar(k -> _pascal_pmf(n, p, k), 0, nothing, x),
        d.n,
        d.p,
        x,
    )
end

function quantile(d::Pascal, u)
    return broadcast(
        (n, p, u) -> _discrete_quantile_scalar(k -> _pascal_pmf(n, p, k), 0, nothing, u),
        d.n,
        d.p,
        u,
    )
end

function rand(rng::AbstractRNG, d::Pascal)
    shape = broadcast_shape(d.n, d.p)
    if isempty(shape)
        return quantile(d, rand(rng))
    end
    u = rand(rng, shape)
    return quantile(d, u)
end

function rand(rng::AbstractRNG, d::Pascal, n::Integer)
    shape = broadcast_shape(d.n, d.p)
    if isempty(shape)
        u = rand(rng, n)
        return quantile(d, u)
    end
    u = rand(rng, (n, shape...))
    return quantile(d, u)
end

mean(d::Pascal) = @. d.n * (1 - d.p) / d.p
var(d::Pascal) = @. d.n * (1 - d.p) / (d.p ^ 2)
