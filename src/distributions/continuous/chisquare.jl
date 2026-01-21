import SpecialFunctions: gamma

struct ChiSquare{T} <: ContinuousDistribution
    n::T
end

support(::ChiSquare) = (0.0, Inf)

function _chisquare_pdf(n::Real, x::Real)
    if x <= 0
        return 0.0
    end
    return x ^ (n / 2 - 1) * exp(-x / 2) / (2 ^ (n / 2) * gamma(n / 2))
end

pdf(d::ChiSquare, x) = broadcast((n, x) -> _chisquare_pdf(n, x), d.n, x)

function cdf(d::ChiSquare, x)
    return broadcast((n, x) -> _continuous_cdf_scalar(t -> _chisquare_pdf(n, t), 0.0, Inf, x), d.n, x)
end

function quantile(d::ChiSquare, u)
    _warn_numeric("quantile", "ChiSquare")
    return broadcast((n, u) -> _continuous_quantile_scalar(t -> _chisquare_pdf(n, t), 0.0, Inf, u), d.n, u)
end

rand(rng::AbstractRNG, d::ChiSquare) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::ChiSquare, n::Integer) = quantile(d, rand(rng, n))

mean(d::ChiSquare) = d.n
var(d::ChiSquare) = @. 2 * d.n
