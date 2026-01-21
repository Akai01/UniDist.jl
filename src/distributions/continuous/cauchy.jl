struct Cauchy{Ta,Talpha} <: ContinuousDistribution
    a::Ta
    alpha::Talpha
end

support(::Cauchy) = (-Inf, Inf)

function _cauchy_pdf(a::Real, alpha::Real, x::Real)
    return 1 / (alpha * pi * (1 + ((x - a) / alpha) ^ 2))
end

pdf(d::Cauchy, x) = broadcast((a, alpha, x) -> _cauchy_pdf(a, alpha, x), d.a, d.alpha, x)

cdf(d::Cauchy, x) = broadcast((a, alpha, x) -> 0.5 + atan((x - a) / alpha) / pi, d.a, d.alpha, x)

quantile(d::Cauchy, u) = broadcast(
    (a, alpha, u) -> u <= 0 ? -Inf : (u >= 1 ? Inf : a + alpha * tan(pi * (u - 0.5))),
    d.a,
    d.alpha,
    u,
)

rand(rng::AbstractRNG, d::Cauchy) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::Cauchy, n::Integer) = quantile(d, rand(rng, n))
