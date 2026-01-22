import SpecialFunctions: beta

struct Beta{Tbeta,Tgamma} <: ContinuousDistribution
    beta::Tbeta
    gamma::Tgamma
end

support(::Beta) = (0.0, 1.0)

function _beta_pdf(beta_param::Real, gamma_param::Real, x::Real)
    if x <= 0 || x >= 1
        return 0.0
    end
    return x ^ (beta_param - 1) * (1 - x) ^ (gamma_param - 1) / beta(beta_param, gamma_param)
end

pdf(d::Beta, x) = broadcast((b, g, x) -> _beta_pdf(b, g, x), d.beta, d.gamma, x)

function cdf(d::Beta, x)
    _info_numeric("cdf", "Beta")
    return broadcast(
        (b, g, x) -> _continuous_cdf_scalar(t -> _beta_pdf(b, g, t), 0.0, 1.0, x),
        d.beta,
        d.gamma,
        x,
    )
end

function quantile(d::Beta, u)
    _info_numeric("quantile", "Beta")
    return broadcast(
        (b, g, u) -> _continuous_quantile_scalar(t -> _beta_pdf(b, g, t), 0.0, 1.0, u),
        d.beta,
        d.gamma,
        u,
    )
end

rand(rng::AbstractRNG, d::Beta) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::Beta, n::Integer) = quantile(d, rand(rng, n))

mean(d::Beta) = @. d.beta / (d.beta + d.gamma)
var(d::Beta) = @. (d.beta * d.gamma) / ((d.beta + d.gamma) ^ 2 * (d.beta + d.gamma + 1))
