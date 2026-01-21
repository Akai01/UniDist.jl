import SpecialFunctions: loggamma

struct Poisson{T} <: DiscreteDistribution
    mu::T
end

support(::Poisson) = nothing

function _poisson_pmf(mu::Real, x::Real)
    if !is_integer_value(x) || x < 0
        return 0.0
    end
    k = Int(x)
    logpmf = -mu + k * log(mu) - loggamma(k + 1)
    return exp(logpmf)
end

pdf(d::Poisson, x) = broadcast((mu, x) -> _poisson_pmf(mu, x), d.mu, x)

function cdf(d::Poisson, x)
    return broadcast(
        (mu, x) -> _discrete_cdf_scalar(k -> _poisson_pmf(mu, k), 0, nothing, x),
        d.mu,
        x,
    )
end

function quantile(d::Poisson, u)
    return broadcast(
        (mu, u) -> _discrete_quantile_scalar(k -> _poisson_pmf(mu, k), 0, nothing, u),
        d.mu,
        u,
    )
end

function rand(rng::AbstractRNG, d::Poisson)
    shape = broadcast_shape(d.mu)
    if isempty(shape)
        return quantile(d, rand(rng))
    end
    u = rand(rng, shape)
    return quantile(d, u)
end

function rand(rng::AbstractRNG, d::Poisson, n::Integer)
    shape = broadcast_shape(d.mu)
    if isempty(shape)
        u = rand(rng, n)
        return quantile(d, u)
    end
    u = rand(rng, (n, shape...))
    return quantile(d, u)
end

mean(d::Poisson) = d.mu
var(d::Poisson) = d.mu
