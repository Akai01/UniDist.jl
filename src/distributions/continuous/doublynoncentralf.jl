import SpecialFunctions: loggamma

struct DoublyNoncentralF{Tn1,Tn2,Tdelta,Tgamma} <: ContinuousDistribution
    n1::Tn1
    n2::Tn2
    delta::Tdelta
    gamma_param::Tgamma
end

support(::DoublyNoncentralF) = (0.0, Inf)

function _doubly_noncentral_f_pdf(n1::Real, n2::Real, delta::Real, gamma_param::Real, x::Real)
    if x <= 0
        return 0.0
    end
    s = 0.0
    terms = min(_SERIES_TERMS, 30)
    for j in 0:terms
        for k in 0:terms
            logw = -delta / 2 + j * log(delta / 2) - loggamma(j + 1) - gamma_param / 2 + k * log(gamma_param / 2) - loggamma(k + 1)
            logbeta = loggamma(n1 / 2 + j) + loggamma(n2 / 2 + k) - loggamma((n1 + n2) / 2 + j + k)
            logterm = logw + (n1 / 2 + j) * log(n1) + (n2 / 2 + k) * log(n2) + (n1 / 2 + j - 1) * log(x) -
                ((n1 + n2) / 2 + j + k) * log(n2 + n1 * x) - logbeta
            s += exp(logterm)
        end
    end
    return s
end

pdf(d::DoublyNoncentralF, x) =
    broadcast((n1, n2, delta, gamma_param, x) -> _doubly_noncentral_f_pdf(n1, n2, delta, gamma_param, x), d.n1, d.n2, d.delta, d.gamma_param, x)

function cdf(d::DoublyNoncentralF, x)
    _warn_numeric("cdf", "DoublyNoncentralF")
    return broadcast(
        (n1, n2, delta, gamma_param, x) -> _continuous_cdf_scalar(t -> _doubly_noncentral_f_pdf(n1, n2, delta, gamma_param, t), 0.0, Inf, x),
        d.n1,
        d.n2,
        d.delta,
        d.gamma_param,
        x,
    )
end

function quantile(d::DoublyNoncentralF, u)
    _warn_numeric("quantile", "DoublyNoncentralF")
    return broadcast(
        (n1, n2, delta, gamma_param, u) -> _continuous_quantile_scalar(t -> _doubly_noncentral_f_pdf(n1, n2, delta, gamma_param, t), 0.0, Inf, u),
        d.n1,
        d.n2,
        d.delta,
        d.gamma_param,
        u,
    )
end

rand(rng::AbstractRNG, d::DoublyNoncentralF) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::DoublyNoncentralF, n::Integer) = quantile(d, rand(rng, n))
