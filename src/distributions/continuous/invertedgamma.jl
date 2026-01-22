import SpecialFunctions: gamma

struct InvertedGamma{Talpha,Tbeta} <: ContinuousDistribution
    alpha::Talpha
    beta::Tbeta
end

support(::InvertedGamma) = (0.0, Inf)

function _inverted_gamma_pdf(alpha::Real, beta::Real, x::Real)
    if x <= 0
        return 0.0
    end
    return x ^ (-(alpha + 1)) * exp(-1 / (beta * x)) / (gamma(alpha) * beta ^ alpha)
end

pdf(d::InvertedGamma, x) = broadcast((alpha, beta, x) -> _inverted_gamma_pdf(alpha, beta, x), d.alpha, d.beta, x)

function cdf(d::InvertedGamma, x)
    _info_numeric("cdf", "InvertedGamma")
    return broadcast(
        (alpha, beta, x) -> _continuous_cdf_scalar(t -> _inverted_gamma_pdf(alpha, beta, t), 0.0, Inf, x),
        d.alpha,
        d.beta,
        x,
    )
end

function quantile(d::InvertedGamma, u)
    _info_numeric("quantile", "InvertedGamma")
    return broadcast(
        (alpha, beta, u) -> _continuous_quantile_scalar(t -> _inverted_gamma_pdf(alpha, beta, t), 0.0, Inf, u),
        d.alpha,
        d.beta,
        u,
    )
end

rand(rng::AbstractRNG, d::InvertedGamma) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::InvertedGamma, n::Integer) = quantile(d, rand(rng, n))

mean(d::InvertedGamma) = @. d.alpha > 1 ? 1 / (d.beta * (d.alpha - 1)) : Inf
var(d::InvertedGamma) = @. d.alpha > 2 ? 1 / (d.beta ^ 2 * (d.alpha - 1) ^ 2 * (d.alpha - 2)) : Inf
