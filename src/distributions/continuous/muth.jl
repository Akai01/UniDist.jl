struct Muth{T} <: ContinuousDistribution
    kappa::T
end

support(::Muth) = (0.0, Inf)

function _muth_survival(kappa::Real, x::Real)
    return exp(1 / kappa - exp(kappa * x) / kappa + kappa * x)
end

function _muth_pdf(kappa::Real, x::Real)
    if x <= 0
        return 0.0
    end
    s = _muth_survival(kappa, x)
    return (exp(kappa * x) - kappa) * s
end

pdf(d::Muth, x) = broadcast((kappa, x) -> _muth_pdf(kappa, x), d.kappa, x)

cdf(d::Muth, x) = broadcast((kappa, x) -> x <= 0 ? 0.0 : 1 - _muth_survival(kappa, x), d.kappa, x)

function quantile(d::Muth, u)
    _warn_numeric("quantile", "Muth")
    return broadcast(
        (kappa, u) -> _continuous_quantile_scalar(t -> _muth_pdf(kappa, t), 0.0, Inf, u),
        d.kappa,
        u,
    )
end

rand(rng::AbstractRNG, d::Muth) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::Muth, n::Integer) = quantile(d, rand(rng, n))
