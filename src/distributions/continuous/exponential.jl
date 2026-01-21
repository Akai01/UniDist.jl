struct Exponential{T} <: ContinuousDistribution
    alpha::T
end

support(::Exponential) = (0.0, Inf)

function _exponential_pdf(alpha::Real, x::Real)
    if x <= 0
        return 0.0
    end
    return (1 / alpha) * exp(-x / alpha)
end

pdf(d::Exponential, x) = broadcast((alpha, x) -> _exponential_pdf(alpha, x), d.alpha, x)

cdf(d::Exponential, x) = broadcast((alpha, x) -> x <= 0 ? 0.0 : 1 - exp(-x / alpha), d.alpha, x)

quantile(d::Exponential, u) = broadcast((alpha, u) -> u <= 0 ? 0.0 : (u >= 1 ? Inf : -alpha * log(1 - u)), d.alpha, u)

rand(rng::AbstractRNG, d::Exponential) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::Exponential, n::Integer) = quantile(d, rand(rng, n))

mean(d::Exponential) = d.alpha
var(d::Exponential) = @. d.alpha ^ 2
