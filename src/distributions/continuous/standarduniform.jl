struct StandardUniform <: ContinuousDistribution end

support(::StandardUniform) = (0.0, 1.0)

pdf(::StandardUniform, x) = @. ifelse((x > 0) & (x < 1), 1.0, 0.0)

cdf(::StandardUniform, x) = @. ifelse(x <= 0, 0.0, ifelse(x >= 1, 1.0, x))

quantile(::StandardUniform, u) = @. ifelse(u <= 0, 0.0, ifelse(u >= 1, 1.0, u))

rand(rng::AbstractRNG, d::StandardUniform) = rand(rng)
rand(rng::AbstractRNG, d::StandardUniform, n::Integer) = rand(rng, n)

mean(::StandardUniform) = 0.5
var(::StandardUniform) = 1 / 12
