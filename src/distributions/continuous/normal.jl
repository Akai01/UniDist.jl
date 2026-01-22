import SpecialFunctions: erf

struct Normal{Tmu,Tsigma} <: ContinuousDistribution
    mu::Tmu
    sigma::Tsigma
end

support(::Normal) = (-Inf, Inf)

function _normal_pdf(mu::Real, sigma::Real, x::Real)
    return exp(-0.5 * ((x - mu) / sigma) ^ 2) / (sigma * sqrt(2 * pi))
end

pdf(d::Normal, x) = broadcast((mu, sigma, x) -> _normal_pdf(mu, sigma, x), d.mu, d.sigma, x)

cdf(d::Normal, x) = broadcast((mu, sigma, x) -> 0.5 * (1 + erf((x - mu) / (sigma * sqrt(2)))), d.mu, d.sigma, x)

function quantile(d::Normal, u)
    _info_numeric("quantile", "Normal")
    return broadcast(
        (mu, sigma, u) -> _continuous_quantile_scalar(t -> _normal_pdf(mu, sigma, t), -Inf, Inf, u),
        d.mu,
        d.sigma,
        u,
    )
end

rand(rng::AbstractRNG, d::Normal) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::Normal, n::Integer) = quantile(d, rand(rng, n))

mean(d::Normal) = d.mu
var(d::Normal) = @. d.sigma ^ 2
