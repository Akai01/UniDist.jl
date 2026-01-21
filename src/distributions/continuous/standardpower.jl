struct StandardPower{T} <: ContinuousDistribution
    beta::T
end

support(::StandardPower) = (0.0, 1.0)

function _standard_power_pdf(beta::Real, x::Real)
    if x <= 0 || x >= 1
        return 0.0
    end
    return beta * x ^ (beta - 1)
end

pdf(d::StandardPower, x) = broadcast((beta, x) -> _standard_power_pdf(beta, x), d.beta, x)

cdf(d::StandardPower, x) = broadcast((beta, x) -> x <= 0 ? 0.0 : (x >= 1 ? 1.0 : x ^ beta), d.beta, x)

quantile(d::StandardPower, u) = broadcast((beta, u) -> u <= 0 ? 0.0 : (u >= 1 ? 1.0 : u ^ (1 / beta)), d.beta, u)

rand(rng::AbstractRNG, d::StandardPower) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::StandardPower, n::Integer) = quantile(d, rand(rng, n))

mean(d::StandardPower) = @. d.beta / (d.beta + 1)
var(d::StandardPower) = @. d.beta / ((d.beta + 1) ^ 2 * (d.beta + 2))
