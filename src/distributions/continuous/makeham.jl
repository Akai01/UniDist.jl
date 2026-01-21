struct Makeham{Tdelta,Tkappa,Tgamma} <: ContinuousDistribution
    delta::Tdelta
    kappa::Tkappa
    gamma_param::Tgamma
end

support(::Makeham) = (0.0, Inf)

function _makeham_pdf(delta::Real, kappa::Real, gamma_param::Real, x::Real)
    if x <= 0
        return 0.0
    end
    return (gamma_param + delta * (kappa ^ x)) * exp(-gamma_param * x - delta * (kappa ^ x - 1) / log(kappa))
end

pdf(d::Makeham, x) =
    broadcast((delta, kappa, gamma_param, x) -> _makeham_pdf(delta, kappa, gamma_param, x), d.delta, d.kappa, d.gamma_param, x)

cdf(d::Makeham, x) =
    broadcast((delta, kappa, gamma_param, x) -> x <= 0 ? 0.0 : 1 - exp(-gamma_param * x - delta * (kappa ^ x - 1) / log(kappa)), d.delta, d.kappa, d.gamma_param, x)

function quantile(d::Makeham, u)
    _warn_numeric("quantile", "Makeham")
    return broadcast(
        (delta, kappa, gamma_param, u) -> _continuous_quantile_scalar(t -> _makeham_pdf(delta, kappa, gamma_param, t), 0.0, Inf, u),
        d.delta,
        d.kappa,
        d.gamma_param,
        u,
    )
end

rand(rng::AbstractRNG, d::Makeham) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::Makeham, n::Integer) = quantile(d, rand(rng, n))
