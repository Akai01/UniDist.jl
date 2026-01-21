struct Lomax{Tlambda,Tkappa} <: ContinuousDistribution
    lambda::Tlambda
    kappa::Tkappa
end

support(::Lomax) = (0.0, Inf)

function _lomax_pdf(lambda::Real, kappa::Real, x::Real)
    if x <= 0
        return 0.0
    end
    return (lambda * kappa) / (1 + lambda * x) ^ (kappa + 1)
end

pdf(d::Lomax, x) = broadcast((lambda, kappa, x) -> _lomax_pdf(lambda, kappa, x), d.lambda, d.kappa, x)

cdf(d::Lomax, x) = broadcast((lambda, kappa, x) -> x <= 0 ? 0.0 : 1 - (1 + lambda * x) ^ (-kappa), d.lambda, d.kappa, x)

quantile(d::Lomax, u) = broadcast((lambda, kappa, u) -> u <= 0 ? 0.0 : (u >= 1 ? Inf : ( (1 - u) ^ (-1 / kappa) - 1 ) / lambda), d.lambda, d.kappa, u)

rand(rng::AbstractRNG, d::Lomax) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::Lomax, n::Integer) = quantile(d, rand(rng, n))

mean(d::Lomax) = @. d.kappa > 1 ? 1 / (d.lambda * (d.kappa - 1)) : Inf
var(d::Lomax) = @. d.kappa > 2 ? d.kappa / (d.lambda ^ 2 * (d.kappa - 1) ^ 2 * (d.kappa - 2)) : Inf
