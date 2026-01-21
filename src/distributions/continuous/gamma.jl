import SpecialFunctions: gamma

struct Gamma{Talpha,Tbeta} <: ContinuousDistribution
    alpha::Talpha
    beta::Tbeta
end

support(::Gamma) = (0.0, Inf)

function _gamma_pdf(alpha::Real, beta::Real, x::Real)
    if x <= 0
        return 0.0
    end
    return x ^ (beta - 1) * exp(-x / alpha) / (alpha ^ beta * gamma(beta))
end

pdf(d::Gamma, x) = broadcast((alpha, beta, x) -> _gamma_pdf(alpha, beta, x), d.alpha, d.beta, x)

function cdf(d::Gamma, x)
    return broadcast(
        (alpha, beta, x) -> _continuous_cdf_scalar(t -> _gamma_pdf(alpha, beta, t), 0.0, Inf, x),
        d.alpha,
        d.beta,
        x,
    )
end

function quantile(d::Gamma, u)
    _warn_numeric("quantile", "Gamma")
    return broadcast(
        (alpha, beta, u) -> _continuous_quantile_scalar(t -> _gamma_pdf(alpha, beta, t), 0.0, Inf, u),
        d.alpha,
        d.beta,
        u,
    )
end

rand(rng::AbstractRNG, d::Gamma) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::Gamma, n::Integer) = quantile(d, rand(rng, n))

mean(d::Gamma) = @. d.alpha * d.beta
var(d::Gamma) = @. d.alpha ^ 2 * d.beta
