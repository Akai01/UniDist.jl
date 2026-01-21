struct Arcsin <: ContinuousDistribution end

support(::Arcsin) = (0.0, 1.0)

function pdf(::Arcsin, x)
    return @. ifelse((x > 0) & (x < 1), 1 / (pi * sqrt(x * (1 - x))), 0.0)
end

function cdf(::Arcsin, x)
    return @. ifelse(x <= 0, 0.0, ifelse(x >= 1, 1.0, (pi + 2 * asin(2 * x - 1)) / (2 * pi)))
end

function quantile(::Arcsin, u)
    return @. ifelse(u <= 0, 0.0, ifelse(u >= 1, 1.0, 0.5 - 0.5 * cos(pi * u)))
end

rand(rng::AbstractRNG, d::Arcsin) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::Arcsin, n::Integer) = quantile(d, rand(rng, n))

mean(::Arcsin) = 0.5
var(::Arcsin) = 0.125
