struct Pareto{Tlambda,Tkappa} <: ContinuousDistribution
    lambda::Tlambda
    kappa::Tkappa
end

support(d::Pareto) = (d.lambda, Inf)

function _pareto_pdf(lambda::Real, kappa::Real, x::Real)
    if x <= lambda
        return 0.0
    end
    return kappa * lambda ^ kappa / x ^ (kappa + 1)
end

pdf(d::Pareto, x) = broadcast((lambda, kappa, x) -> _pareto_pdf(lambda, kappa, x), d.lambda, d.kappa, x)

cdf(d::Pareto, x) = broadcast((lambda, kappa, x) -> x <= lambda ? 0.0 : 1 - (lambda / x) ^ kappa, d.lambda, d.kappa, x)

quantile(d::Pareto, u) = broadcast((lambda, kappa, u) -> u <= 0 ? lambda : (u >= 1 ? Inf : lambda / (1 - u) ^ (1 / kappa)), d.lambda, d.kappa, u)

rand(rng::AbstractRNG, d::Pareto) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::Pareto, n::Integer) = quantile(d, rand(rng, n))

mean(d::Pareto) = @. d.kappa > 1 ? d.kappa * d.lambda / (d.kappa - 1) : Inf
var(d::Pareto) = @. d.kappa > 2 ? (d.kappa * d.lambda ^ 2) / ((d.kappa - 1) ^ 2 * (d.kappa - 2)) : Inf
