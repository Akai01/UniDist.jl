struct GeneralizedPareto{Tdelta,Tkappa,Tgamma} <: ContinuousDistribution
    delta::Tdelta
    kappa::Tkappa
    gamma_param::Tgamma
end

support(::GeneralizedPareto) = (0.0, Inf)

function _generalized_pareto_pdf(delta::Real, kappa::Real, gamma_param::Real, x::Real)
    if x <= 0
        return 0.0
    end
    return (gamma_param + kappa / (x + delta)) * (1 + x / delta) ^ (-kappa) * exp(-gamma_param * x)
end

pdf(d::GeneralizedPareto, x) =
    broadcast((delta, kappa, gamma_param, x) -> _generalized_pareto_pdf(delta, kappa, gamma_param, x), d.delta, d.kappa, d.gamma_param, x)

cdf(d::GeneralizedPareto, x) =
    broadcast((delta, kappa, gamma_param, x) -> x <= 0 ? 0.0 : 1 - exp(-gamma_param * x) * (1 + x / delta) ^ (-kappa), d.delta, d.kappa, d.gamma_param, x)

function quantile(d::GeneralizedPareto, u)
    _info_numeric("quantile", "GeneralizedPareto")
    return broadcast(
        (delta, kappa, gamma_param, u) -> _continuous_quantile_scalar(
            t -> _generalized_pareto_pdf(delta, kappa, gamma_param, t),
            0.0,
            Inf,
            u,
        ),
        d.delta,
        d.kappa,
        d.gamma_param,
        u,
    )
end

rand(rng::AbstractRNG, d::GeneralizedPareto) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::GeneralizedPareto, n::Integer) = quantile(d, rand(rng, n))
