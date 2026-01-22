import SpecialFunctions: loggamma, gamma

struct NoncentralT{Tn,Tdelta} <: ContinuousDistribution
    n::Tn
    delta::Tdelta
end

support(::NoncentralT) = (-Inf, Inf)

function _noncentral_t_pdf(n::Real, delta::Real, x::Real)
    coef = (n ^ (n / 2)) * exp(-delta ^ 2 / 2) / (sqrt(pi) * gamma(n / 2) * (n + x ^ 2) ^ ((n + 1) / 2))
    s = 0.0
    base = x * delta * sqrt(2) / sqrt(n + x ^ 2)
    for i in 0:_SERIES_TERMS
        logterm = loggamma((n + i + 1) / 2) - loggamma(i + 1) + i * log(abs(base))
        term = exp(logterm)
        if base < 0 && isodd(i)
            term = -term
        end
        s += term
    end
    return coef * s
end

pdf(d::NoncentralT, x) = broadcast((n, delta, x) -> _noncentral_t_pdf(n, delta, x), d.n, d.delta, x)

function cdf(d::NoncentralT, x)
    _info_numeric("cdf", "NoncentralT")
    return broadcast(
        (n, delta, x) -> _continuous_cdf_scalar(t -> _noncentral_t_pdf(n, delta, t), -Inf, Inf, x),
        d.n,
        d.delta,
        x,
    )
end

function quantile(d::NoncentralT, u)
    _info_numeric("quantile", "NoncentralT")
    return broadcast(
        (n, delta, u) -> _continuous_quantile_scalar(t -> _noncentral_t_pdf(n, delta, t), -Inf, Inf, u),
        d.n,
        d.delta,
        u,
    )
end

rand(rng::AbstractRNG, d::NoncentralT) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::NoncentralT, n::Integer) = quantile(d, rand(rng, n))

mean(d::NoncentralT) = @. d.n > 1 ? d.delta * sqrt(d.n / 2) * gamma((d.n - 1) / 2) / gamma(d.n / 2) : NaN

function var(d::NoncentralT)
    return @. d.n > 2 ?
        (d.n * (1 + d.delta ^ 2)) / (d.n - 2) -
        (d.delta ^ 2) * (d.n / 2) * (gamma((d.n - 1) / 2) / gamma(d.n / 2)) ^ 2 :
        NaN
end
