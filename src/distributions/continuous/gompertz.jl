struct Gompertz{Tdelta,Tkappa} <: ContinuousDistribution
    delta::Tdelta
    kappa::Tkappa
end

support(::Gompertz) = (0.0, Inf)

function _gompertz_pdf(delta::Real, kappa::Real, x::Real)
    if x <= 0
        return 0.0
    end
    return delta * (kappa ^ x) * exp(-delta * (kappa ^ x - 1) / log(kappa))
end

pdf(d::Gompertz, x) = broadcast((delta, kappa, x) -> _gompertz_pdf(delta, kappa, x), d.delta, d.kappa, x)

cdf(d::Gompertz, x) = broadcast((delta, kappa, x) -> x <= 0 ? 0.0 : 1 - exp(-delta * (kappa ^ x - 1) / log(kappa)), d.delta, d.kappa, x)

function quantile(d::Gompertz, u)
    _warn_numeric("quantile", "Gompertz")
    return broadcast(
        (delta, kappa, u) -> _continuous_quantile_scalar(t -> _gompertz_pdf(delta, kappa, t), 0.0, Inf, u),
        d.delta,
        d.kappa,
        u,
    )
end

rand(rng::AbstractRNG, d::Gompertz) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::Gompertz, n::Integer) = quantile(d, rand(rng, n))
