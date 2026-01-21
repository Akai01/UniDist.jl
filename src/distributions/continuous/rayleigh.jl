struct Rayleigh{T} <: ContinuousDistribution
    alpha::T
end

support(::Rayleigh) = (0.0, Inf)

function _rayleigh_pdf(alpha::Real, x::Real)
    if x <= 0
        return 0.0
    end
    return (2 * x / alpha) * exp(-(x ^ 2) / alpha)
end

pdf(d::Rayleigh, x) = broadcast((alpha, x) -> _rayleigh_pdf(alpha, x), d.alpha, x)

cdf(d::Rayleigh, x) = broadcast((alpha, x) -> x <= 0 ? 0.0 : 1 - exp(-(x ^ 2) / alpha), d.alpha, x)

quantile(d::Rayleigh, u) = broadcast((alpha, u) -> u <= 0 ? 0.0 : (u >= 1 ? Inf : sqrt(-alpha * log(1 - u))), d.alpha, u)

rand(rng::AbstractRNG, d::Rayleigh) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::Rayleigh, n::Integer) = quantile(d, rand(rng, n))

mean(d::Rayleigh) = @. 0.5 * sqrt(pi * d.alpha)
var(d::Rayleigh) = @. (4 - pi) * d.alpha / 4
