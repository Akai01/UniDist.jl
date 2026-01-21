import SpecialFunctions: loggamma

struct NoncentralF{Tn1,Tn2,Tdelta} <: ContinuousDistribution
    n1::Tn1
    n2::Tn2
    delta::Tdelta
end

support(::NoncentralF) = (0.0, Inf)

function _noncentral_f_pdf(n1::Real, n2::Real, delta::Real, x::Real)
    if x <= 0
        return 0.0
    end
    s = 0.0
    for i in 0:_SERIES_TERMS
        logcoef = loggamma((2 * i + n1 + n2) / 2) - loggamma(n2 / 2) - loggamma(i + n1 / 2) - loggamma(i + 1)
        logterm = logcoef + (2 * i + n1) / 2 * log(n1 / n2) + ((2 * i + n1 - 2) / 2) * log(x) -
            ((2 * i + n1 + n2) / 2) * log(1 + (n1 / n2) * x) - delta / 2 + i * log(delta / 2)
        s += exp(logterm)
    end
    return s
end

pdf(d::NoncentralF, x) = broadcast((n1, n2, delta, x) -> _noncentral_f_pdf(n1, n2, delta, x), d.n1, d.n2, d.delta, x)

function cdf(d::NoncentralF, x)
    _warn_numeric("cdf", "NoncentralF")
    return broadcast(
        (n1, n2, delta, x) -> _continuous_cdf_scalar(t -> _noncentral_f_pdf(n1, n2, delta, t), 0.0, Inf, x),
        d.n1,
        d.n2,
        d.delta,
        x,
    )
end

function quantile(d::NoncentralF, u)
    _warn_numeric("quantile", "NoncentralF")
    return broadcast(
        (n1, n2, delta, u) -> _continuous_quantile_scalar(t -> _noncentral_f_pdf(n1, n2, delta, t), 0.0, Inf, u),
        d.n1,
        d.n2,
        d.delta,
        u,
    )
end

rand(rng::AbstractRNG, d::NoncentralF) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::NoncentralF, n::Integer) = quantile(d, rand(rng, n))

mean(d::NoncentralF) = @. d.n2 > 2 ? d.n2 * (d.n1 + d.delta) / (d.n1 * (d.n2 - 2)) : Inf

function var(d::NoncentralF)
    return @. d.n2 > 4 ?
        (2 * (d.n1 + d.delta) ^ 2 + (d.n1 + 2 * d.delta) * (d.n2 - 2)) * d.n2 ^ 2 / (d.n1 ^ 2 * (d.n2 - 2) ^ 2 * (d.n2 - 4)) :
        Inf
end
