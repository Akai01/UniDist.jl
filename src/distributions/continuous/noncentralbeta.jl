import SpecialFunctions: loggamma

struct NoncentralBeta{Tbeta,Tgamma,Tdelta} <: ContinuousDistribution
    beta::Tbeta
    gamma::Tgamma
    delta::Tdelta
end

support(::NoncentralBeta) = (0.0, 1.0)

function _noncentral_beta_pdf(beta_param::Real, gamma_param::Real, delta::Real, x::Real)
    if x <= 0 || x >= 1
        return 0.0
    end
    s = 0.0
    for i in 0:_SERIES_TERMS
        logcoef = loggamma(i + beta_param + gamma_param) - loggamma(gamma_param) - loggamma(i + beta_param) - loggamma(i + 1)
        logterm = logcoef - delta / 2 + i * log(delta / 2) + (i + beta_param - 1) * log(x) + (gamma_param - 1) * log(1 - x)
        s += exp(logterm)
    end
    return s
end

pdf(d::NoncentralBeta, x) = broadcast((beta_param, gamma_param, delta, x) -> _noncentral_beta_pdf(beta_param, gamma_param, delta, x), d.beta, d.gamma, d.delta, x)

function cdf(d::NoncentralBeta, x)
    _info_numeric("cdf", "NoncentralBeta")
    return broadcast(
        (beta_param, gamma_param, delta, x) -> _continuous_cdf_scalar(t -> _noncentral_beta_pdf(beta_param, gamma_param, delta, t), 0.0, 1.0, x),
        d.beta,
        d.gamma,
        d.delta,
        x,
    )
end

function quantile(d::NoncentralBeta, u)
    _info_numeric("quantile", "NoncentralBeta")
    return broadcast(
        (beta_param, gamma_param, delta, u) -> _continuous_quantile_scalar(t -> _noncentral_beta_pdf(beta_param, gamma_param, delta, t), 0.0, 1.0, u),
        d.beta,
        d.gamma,
        d.delta,
        u,
    )
end

rand(rng::AbstractRNG, d::NoncentralBeta) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::NoncentralBeta, n::Integer) = quantile(d, rand(rng, n))
