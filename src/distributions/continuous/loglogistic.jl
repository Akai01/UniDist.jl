struct LogLogistic{Tlambda,Tkappa} <: ContinuousDistribution
    lambda::Tlambda
    kappa::Tkappa
end

support(::LogLogistic) = (0.0, Inf)

function _loglogistic_pdf(lambda::Real, kappa::Real, x::Real)
    if x <= 0
        return 0.0
    end
    num = lambda ^ kappa * x ^ (kappa - 1)
    den = (1 + (lambda * x) ^ kappa) ^ 2
    return kappa * num / den
end

pdf(d::LogLogistic, x) = broadcast((lambda, kappa, x) -> _loglogistic_pdf(lambda, kappa, x), d.lambda, d.kappa, x)

cdf(d::LogLogistic, x) = broadcast((lambda, kappa, x) -> x <= 0 ? 0.0 : (lambda * x) ^ kappa / (1 + (lambda * x) ^ kappa), d.lambda, d.kappa, x)

quantile(d::LogLogistic, u) = broadcast((lambda, kappa, u) -> u <= 0 ? 0.0 : (u >= 1 ? Inf : (u / (1 - u)) ^ (1 / kappa) / lambda), d.lambda, d.kappa, u)

rand(rng::AbstractRNG, d::LogLogistic) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::LogLogistic, n::Integer) = quantile(d, rand(rng, n))
