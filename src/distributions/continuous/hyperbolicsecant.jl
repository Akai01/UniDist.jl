struct HyperbolicSecant <: ContinuousDistribution end

support(::HyperbolicSecant) = (-Inf, Inf)

sech(x::Real) = 1 / cosh(x)

pdf(::HyperbolicSecant, x) = @. sech(pi * x)

cdf(::HyperbolicSecant, x) = @. (pi + 2 * atan(sinh(pi * x))) / (2 * pi)

quantile(::HyperbolicSecant, u) = @. asinh(tan(pi * (u - 0.5))) / pi

rand(rng::AbstractRNG, d::HyperbolicSecant) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::HyperbolicSecant, n::Integer) = quantile(d, rand(rng, n))
