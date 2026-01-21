struct ExtremeValue{Talpha,Tbeta} <: ContinuousDistribution
    alpha::Talpha
    beta::Tbeta
end

support(::ExtremeValue) = (-Inf, Inf)

# Paper parameterization: f(x) = (β/α) exp(βx - e^(βx)/α), -∞ < x < ∞
function _extreme_value_pdf(alpha::Real, beta::Real, x::Real)
    z = exp(beta * x)
    return (beta / alpha) * z * exp(-z / alpha)
end

pdf(d::ExtremeValue, x) = broadcast((alpha, beta, x) -> _extreme_value_pdf(alpha, beta, x), d.alpha, d.beta, x)

cdf(d::ExtremeValue, x) = broadcast((alpha, beta, x) -> 1 - exp(-exp(beta * x) / alpha), d.alpha, d.beta, x)

quantile(d::ExtremeValue, u) = broadcast((alpha, beta, u) -> u <= 0 ? -Inf : (u >= 1 ? Inf : log(-alpha * log(1 - u)) / beta), d.alpha, d.beta, u)

rand(rng::AbstractRNG, d::ExtremeValue) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::ExtremeValue, n::Integer) = quantile(d, rand(rng, n))
