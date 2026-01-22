import SpecialFunctions: loggamma

struct DoublyNoncentralT{Tn,Tdelta,Tgamma} <: ContinuousDistribution
    n::Tn
    delta::Tdelta
    gamma_param::Tgamma
end

support(::DoublyNoncentralT) = (-Inf, Inf)

function _doubly_noncentral_t_pdf(n::Real, delta::Real, gamma_param::Real, x::Real)
    s = 0.0
    terms = min(_SERIES_TERMS, 50)
    for j in 0:terms
        logw = -gamma_param / 2 + j * log(gamma_param / 2) - loggamma(j + 1)
        w = exp(logw)
        s += w * _noncentral_t_pdf(n + 2 * j, delta, x)
    end
    return s
end

pdf(d::DoublyNoncentralT, x) = broadcast((n, delta, gamma_param, x) -> _doubly_noncentral_t_pdf(n, delta, gamma_param, x), d.n, d.delta, d.gamma_param, x)

function cdf(d::DoublyNoncentralT, x)
    _info_numeric("cdf", "DoublyNoncentralT")
    return broadcast(
        (n, delta, gamma_param, x) -> _continuous_cdf_scalar(t -> _doubly_noncentral_t_pdf(n, delta, gamma_param, t), -Inf, Inf, x),
        d.n,
        d.delta,
        d.gamma_param,
        x,
    )
end

function quantile(d::DoublyNoncentralT, u)
    _info_numeric("quantile", "DoublyNoncentralT")
    return broadcast(
        (n, delta, gamma_param, u) -> _continuous_quantile_scalar(t -> _doubly_noncentral_t_pdf(n, delta, gamma_param, t), -Inf, Inf, u),
        d.n,
        d.delta,
        d.gamma_param,
        u,
    )
end

rand(rng::AbstractRNG, d::DoublyNoncentralT) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::DoublyNoncentralT, n::Integer) = quantile(d, rand(rng, n))
