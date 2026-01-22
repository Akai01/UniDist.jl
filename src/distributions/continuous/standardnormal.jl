import SpecialFunctions: erf

struct StandardNormal <: ContinuousDistribution end

support(::StandardNormal) = (-Inf, Inf)

pdf(::StandardNormal, x) = @. exp(-0.5 * x ^ 2) / sqrt(2 * pi)

cdf(::StandardNormal, x) = @. 0.5 * (1 + erf(x / sqrt(2)))

function quantile(::StandardNormal, u)
    _info_numeric("quantile", "StandardNormal")
    return broadcast(u -> _continuous_quantile_scalar(t -> exp(-0.5 * t ^ 2) / sqrt(2 * pi), -Inf, Inf, u), u)
end

rand(rng::AbstractRNG, d::StandardNormal) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::StandardNormal, n::Integer) = quantile(d, rand(rng, n))

mean(::StandardNormal) = 0.0
var(::StandardNormal) = 1.0
