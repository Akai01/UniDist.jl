import SpecialFunctions: loggamma

struct NoncentralChiSquare{Tdelta,Tn} <: ContinuousDistribution
    delta::Tdelta
    n::Tn
end

support(::NoncentralChiSquare) = (0.0, Inf)

function _noncentral_chisquare_pdf(delta::Real, n::Real, x::Real)
    if x <= 0
        return 0.0
    end
    s = 0.0
    for k in 0:_SERIES_TERMS
        logterm = -delta / 2 + k * log(delta / 2) - loggamma(k + 1) - x / 2 +
            (n / 2 + k - 1) * log(x) - (n / 2 + k) * log(2) - loggamma(n / 2 + k)
        s += exp(logterm)
    end
    return s
end

pdf(d::NoncentralChiSquare, x) = broadcast((delta, n, x) -> _noncentral_chisquare_pdf(delta, n, x), d.delta, d.n, x)

function cdf(d::NoncentralChiSquare, x)
    _warn_numeric("cdf", "NoncentralChiSquare")
    return broadcast(
        (delta, n, x) -> _continuous_cdf_scalar(t -> _noncentral_chisquare_pdf(delta, n, t), 0.0, Inf, x),
        d.delta,
        d.n,
        x,
    )
end

function quantile(d::NoncentralChiSquare, u)
    _warn_numeric("quantile", "NoncentralChiSquare")
    return broadcast(
        (delta, n, u) -> _continuous_quantile_scalar(t -> _noncentral_chisquare_pdf(delta, n, t), 0.0, Inf, u),
        d.delta,
        d.n,
        u,
    )
end

rand(rng::AbstractRNG, d::NoncentralChiSquare) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::NoncentralChiSquare, n::Integer) = quantile(d, rand(rng, n))

mean(d::NoncentralChiSquare) = @. d.delta + d.n
var(d::NoncentralChiSquare) = @. 2 * (d.n + 2 * d.delta)
