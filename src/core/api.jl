import Random: AbstractRNG, default_rng, rand
import Statistics: mean, var

function pdf end
function cdf end
function quantile end

# Shorthand API (R-like)
d(dist::AbstractDistribution, x) = pdf(dist, x)
p(dist::AbstractDistribution, x) = cdf(dist, x)
q(dist::AbstractDistribution, u) = quantile(dist, u)

r(dist::AbstractDistribution, n::Integer=1; rng::AbstractRNG=default_rng()) = rand(rng, dist, n)

sf(dist::AbstractDistribution, x) = 1 .- cdf(dist, x)
hazard(dist::AbstractDistribution, x) = pdf(dist, x) ./ sf(dist, x)
cumhaz(dist::AbstractDistribution, x) = .-log.(sf(dist, x))

mean(dist::AbstractDistribution) = throw(MethodError(mean, (dist,)))
var(dist::AbstractDistribution) = throw(MethodError(var, (dist,)))
