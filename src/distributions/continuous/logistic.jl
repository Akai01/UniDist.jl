struct Logistic{Tlambda,Tkappa} <: ContinuousDistribution
    lambda::Tlambda
    kappa::Tkappa
end

support(::Logistic) = (-Inf, Inf)

function _logistic_pdf(lambda::Real, kappa::Real, x::Real)
    num = (lambda ^ kappa) * kappa * exp(kappa * x)
    den = (1 + (lambda * exp(x)) ^ kappa) ^ 2
    return num / den
end

pdf(d::Logistic, x) = broadcast((lambda, kappa, x) -> _logistic_pdf(lambda, kappa, x), d.lambda, d.kappa, x)

cdf(d::Logistic, x) = broadcast((lambda, kappa, x) -> (lambda ^ kappa) * exp(kappa * x) / (1 + (lambda ^ kappa) * exp(kappa * x)), d.lambda, d.kappa, x)

quantile(d::Logistic, u) = broadcast((lambda, kappa, u) -> begin
    if u <= 0
        -Inf
    elseif u >= 1
        Inf
    else
        (log(u / (1 - u)) / kappa) - log(lambda)
    end
end, d.lambda, d.kappa, u)

rand(rng::AbstractRNG, d::Logistic) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::Logistic, n::Integer) = quantile(d, rand(rng, n))
