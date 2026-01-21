struct StandardCauchy <: ContinuousDistribution end

support(::StandardCauchy) = (-Inf, Inf)

pdf(::StandardCauchy, x) = @. 1 / (pi * (1 + x ^ 2))

cdf(::StandardCauchy, x) = @. 0.5 + atan(x) / pi

quantile(::StandardCauchy, u) = @. ifelse(u <= 0, -Inf, ifelse(u >= 1, Inf, tan(pi * (u - 0.5))))

rand(rng::AbstractRNG, d::StandardCauchy) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::StandardCauchy, n::Integer) = quantile(d, rand(rng, n))
