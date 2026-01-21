import SpecialFunctions: gamma

struct Error{Ta,Tb,Tc} <: ContinuousDistribution
    a::Ta
    b::Tb
    c::Tc
end

support(::Error) = (-Inf, Inf)

function _error_pdf(a::Real, b::Real, c::Real, x::Real)
    z = abs(x - a) / b
    return exp(-0.5 * z ^ (2 / c)) / (b * (2 ^ (c / 2 + 1)) * gamma(c / 2 + 1))
end

pdf(d::Error, x) = broadcast((a, b, c, x) -> _error_pdf(a, b, c, x), d.a, d.b, d.c, x)

function cdf(d::Error, x)
    _warn_numeric("cdf", "Error")
    return broadcast(
        (a, b, c, x) -> _continuous_cdf_scalar(t -> _error_pdf(a, b, c, t), -Inf, Inf, x),
        d.a,
        d.b,
        d.c,
        x,
    )
end

function quantile(d::Error, u)
    _warn_numeric("quantile", "Error")
    return broadcast(
        (a, b, c, u) -> _continuous_quantile_scalar(t -> _error_pdf(a, b, c, t), -Inf, Inf, u),
        d.a,
        d.b,
        d.c,
        u,
    )
end

rand(rng::AbstractRNG, d::Error) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::Error, n::Integer) = quantile(d, rand(rng, n))
