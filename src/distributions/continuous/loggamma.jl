import SpecialFunctions: gamma

struct LogGamma{Talpha,Tbeta} <: ContinuousDistribution
    alpha::Talpha
    beta::Tbeta
end

support(::LogGamma) = (-Inf, Inf)

function _loggamma_pdf(alpha::Real, beta::Real, x::Real)
    return exp(beta * x - exp(x) / alpha) / (alpha ^ beta * gamma(beta))
end

pdf(d::LogGamma, x) = broadcast((alpha, beta, x) -> _loggamma_pdf(alpha, beta, x), d.alpha, d.beta, x)

function cdf(d::LogGamma, x)
    _warn_numeric("cdf", "LogGamma")
    return broadcast(
        (alpha, beta, x) -> _continuous_cdf_scalar(t -> _loggamma_pdf(alpha, beta, t), -Inf, Inf, x),
        d.alpha,
        d.beta,
        x,
    )
end

function quantile(d::LogGamma, u)
    _warn_numeric("quantile", "LogGamma")
    return broadcast(
        (alpha, beta, u) -> _continuous_quantile_scalar(t -> _loggamma_pdf(alpha, beta, t), -Inf, Inf, u),
        d.alpha,
        d.beta,
        u,
    )
end

rand(rng::AbstractRNG, d::LogGamma) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::LogGamma, n::Integer) = quantile(d, rand(rng, n))
