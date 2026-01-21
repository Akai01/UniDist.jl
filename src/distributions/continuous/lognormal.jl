import SpecialFunctions: erf, erfinv

struct LogNormal{Talpha,Tbeta} <: ContinuousDistribution
    alpha::Talpha
    beta::Tbeta
end

support(::LogNormal) = (0.0, Inf)

function _lognormal_pdf(alpha::Real, beta::Real, x::Real)
    if x <= 0
        return 0.0
    end
    z = log(x / alpha) / beta
    return exp(-0.5 * z ^ 2) / (x * beta * sqrt(2 * pi))
end

pdf(d::LogNormal, x) = broadcast((alpha, beta, x) -> _lognormal_pdf(alpha, beta, x), d.alpha, d.beta, x)

cdf(d::LogNormal, x) = broadcast((alpha, beta, x) -> x <= 0 ? 0.0 : 0.5 * (1 + erf(log(x / alpha) / (beta * sqrt(2)))), d.alpha, d.beta, x)

quantile(d::LogNormal, u) = broadcast((alpha, beta, u) -> begin
    if u <= 0
        0.0
    elseif u >= 1
        Inf
    else
        alpha * exp(beta * sqrt(2) * erfinv(2 * u - 1))
    end
end, d.alpha, d.beta, u)

rand(rng::AbstractRNG, d::LogNormal) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::LogNormal, n::Integer) = quantile(d, rand(rng, n))

mean(d::LogNormal) = @. d.alpha * exp(d.beta ^ 2 / 2)
var(d::LogNormal) = @. d.alpha ^ 2 * (exp(d.beta ^ 2) - 1) * exp(d.beta ^ 2)
