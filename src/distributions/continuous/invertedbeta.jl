import SpecialFunctions: beta

struct InvertedBeta{Tbeta,Tgamma} <: ContinuousDistribution
    beta::Tbeta
    gamma::Tgamma
end

support(::InvertedBeta) = (0.0, Inf)

function _inverted_beta_pdf(beta_param::Real, gamma_param::Real, x::Real)
    if x <= 0
        return 0.0
    end
    return x ^ (beta_param - 1) * (x + 1) ^ (-beta_param - gamma_param) / beta(beta_param, gamma_param)
end

pdf(d::InvertedBeta, x) = broadcast((b, g, x) -> _inverted_beta_pdf(b, g, x), d.beta, d.gamma, x)

function cdf(d::InvertedBeta, x)
    _warn_numeric("cdf", "InvertedBeta")
    return broadcast(
        (b, g, x) -> _continuous_cdf_scalar(t -> _inverted_beta_pdf(b, g, t), 0.0, Inf, x),
        d.beta,
        d.gamma,
        x,
    )
end

function quantile(d::InvertedBeta, u)
    _warn_numeric("quantile", "InvertedBeta")
    return broadcast(
        (b, g, u) -> _continuous_quantile_scalar(t -> _inverted_beta_pdf(b, g, t), 0.0, Inf, u),
        d.beta,
        d.gamma,
        u,
    )
end

rand(rng::AbstractRNG, d::InvertedBeta) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::InvertedBeta, n::Integer) = quantile(d, rand(rng, n))
