struct ExponentialPower{Tlambda,Tkappa} <: ContinuousDistribution
    lambda::Tlambda
    kappa::Tkappa
end

support(::ExponentialPower) = (0.0, Inf)

function _exponential_power_pdf(lambda::Real, kappa::Real, x::Real)
    if x <= 0
        return 0.0
    end
    t = lambda * x ^ kappa
    return exp(1 - exp(t)) * exp(t) * lambda * kappa * x ^ (kappa - 1)
end

pdf(d::ExponentialPower, x) = broadcast((lambda, kappa, x) -> _exponential_power_pdf(lambda, kappa, x), d.lambda, d.kappa, x)

cdf(d::ExponentialPower, x) = broadcast((lambda, kappa, x) -> x <= 0 ? 0.0 : 1 - exp(1 - exp(lambda * x ^ kappa)), d.lambda, d.kappa, x)

function quantile(d::ExponentialPower, u)
    _info_numeric("quantile", "ExponentialPower")
    return broadcast(
        (lambda, kappa, u) -> _continuous_quantile_scalar(t -> _exponential_power_pdf(lambda, kappa, t), 0.0, Inf, u),
        d.lambda,
        d.kappa,
        u,
    )
end

rand(rng::AbstractRNG, d::ExponentialPower) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::ExponentialPower, n::Integer) = quantile(d, rand(rng, n))
