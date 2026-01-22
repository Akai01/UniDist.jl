import SpecialFunctions: gamma

struct GeneralizedGamma{Talpha,Tbeta,Tgamma} <: ContinuousDistribution
    alpha::Talpha
    beta::Tbeta
    gamma_param::Tgamma
end

support(::GeneralizedGamma) = (0.0, Inf)

function _generalized_gamma_pdf(alpha::Real, beta::Real, gamma_param::Real, x::Real)
    if x <= 0
        return 0.0
    end
    num = gamma_param * x ^ (gamma_param * beta - 1) * exp(-(x / alpha) ^ gamma_param)
    den = alpha ^ (gamma_param * beta) * gamma(beta)
    return num / den
end

pdf(d::GeneralizedGamma, x) = broadcast((alpha, beta, gamma_param, x) -> _generalized_gamma_pdf(alpha, beta, gamma_param, x), d.alpha, d.beta, d.gamma_param, x)

function cdf(d::GeneralizedGamma, x)
    _info_numeric("cdf", "GeneralizedGamma")
    return broadcast(
        (alpha, beta, gamma_param, x) -> _continuous_cdf_scalar(t -> _generalized_gamma_pdf(alpha, beta, gamma_param, t), 0.0, Inf, x),
        d.alpha,
        d.beta,
        d.gamma_param,
        x,
    )
end

function quantile(d::GeneralizedGamma, u)
    _info_numeric("quantile", "GeneralizedGamma")
    return broadcast(
        (alpha, beta, gamma_param, u) -> _continuous_quantile_scalar(t -> _generalized_gamma_pdf(alpha, beta, gamma_param, t), 0.0, Inf, u),
        d.alpha,
        d.beta,
        d.gamma_param,
        u,
    )
end

rand(rng::AbstractRNG, d::GeneralizedGamma) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::GeneralizedGamma, n::Integer) = quantile(d, rand(rng, n))

mean(d::GeneralizedGamma) = @. d.alpha * gamma(d.beta + 1 / d.gamma_param) / gamma(d.beta)
var(d::GeneralizedGamma) = @. d.alpha ^ 2 * (gamma(d.beta + 2 / d.gamma_param) / gamma(d.beta) - (gamma(d.beta + 1 / d.gamma_param) / gamma(d.beta)) ^ 2)
