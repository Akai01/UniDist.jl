struct Power{Ta,Tb} <: ContinuousDistribution
    alpha::Ta
    beta::Tb
end

support(d::Power) = (0.0, d.alpha)

function _power_pdf(alpha::Real, beta::Real, x::Real)
    if x <= 0 || x >= alpha
        return 0.0
    end
    return (beta * x ^ (beta - 1)) / (alpha ^ beta)
end

pdf(d::Power, x) = broadcast((alpha, beta, x) -> _power_pdf(alpha, beta, x), d.alpha, d.beta, x)

function cdf(d::Power, x)
    return broadcast((alpha, beta, x) -> x <= 0 ? 0.0 : (x >= alpha ? 1.0 : (x ^ beta) / (alpha ^ beta)), d.alpha, d.beta, x)
end

function quantile(d::Power, u)
    return broadcast((alpha, beta, u) -> u <= 0 ? 0.0 : (u >= 1 ? alpha : alpha * (u ^ (1 / beta))), d.alpha, d.beta, u)
end

rand(rng::AbstractRNG, d::Power) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::Power, n::Integer) = quantile(d, rand(rng, n))

mean(d::Power) = @. d.alpha * d.beta / (d.beta + 1)
var(d::Power) = @. (d.alpha ^ 2) * d.beta / ((d.beta + 1) ^ 2 * (d.beta + 2))
