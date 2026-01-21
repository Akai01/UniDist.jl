import SpecialFunctions: gamma

struct Weibull{Ta,Tb} <: ContinuousDistribution
    alpha::Ta
    beta::Tb
end

support(::Weibull) = (0.0, Inf)

# Paper parameterization: f(x) = (β/α) x^(β-1) exp[-(1/α)x^β], x > 0
function _weibull_pdf(alpha::Real, beta::Real, x::Real)
    if x <= 0
        return 0.0
    end
    return (beta / alpha) * x ^ (beta - 1) * exp(-(x ^ beta) / alpha)
end

pdf(d::Weibull, x) = broadcast((alpha, beta, x) -> _weibull_pdf(alpha, beta, x), d.alpha, d.beta, x)

cdf(d::Weibull, x) = broadcast((alpha, beta, x) -> x <= 0 ? 0.0 : 1 - exp(-(x ^ beta) / alpha), d.alpha, d.beta, x)

quantile(d::Weibull, u) = broadcast((alpha, beta, u) -> u <= 0 ? 0.0 : (u >= 1 ? Inf : (-alpha * log(1 - u)) ^ (1 / beta)), d.alpha, d.beta, u)

rand(rng::AbstractRNG, d::Weibull) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::Weibull, n::Integer) = quantile(d, rand(rng, n))

mean(d::Weibull) = @. d.alpha ^ (1 / d.beta) * gamma(1 + 1 / d.beta)
var(d::Weibull) = @. d.alpha ^ (2 / d.beta) * (gamma(1 + 2 / d.beta) - gamma(1 + 1 / d.beta) ^ 2)
