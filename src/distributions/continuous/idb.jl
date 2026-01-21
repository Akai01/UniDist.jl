struct IDB{Tdelta,Tkappa,Tgamma} <: ContinuousDistribution
    delta::Tdelta
    kappa::Tkappa
    gamma_param::Tgamma
end

support(::IDB) = (0.0, Inf)

function _idb_pdf(delta::Real, kappa::Real, gamma_param::Real, x::Real)
    if x <= 0
        return 0.0
    end
    num = ((1 + kappa * x) * delta * x + gamma_param) * exp(-delta * x ^ 2 / 2)
    den = (1 + kappa * x) ^ (gamma_param / kappa + 1)
    return num / den
end

pdf(d::IDB, x) = broadcast((delta, kappa, gamma_param, x) -> _idb_pdf(delta, kappa, gamma_param, x), d.delta, d.kappa, d.gamma_param, x)

function cdf(d::IDB, x)
    _warn_numeric("cdf", "IDB")
    return broadcast(
        (delta, kappa, gamma_param, x) -> _continuous_cdf_scalar(t -> _idb_pdf(delta, kappa, gamma_param, t), 0.0, Inf, x),
        d.delta,
        d.kappa,
        d.gamma_param,
        x,
    )
end

function quantile(d::IDB, u)
    _warn_numeric("quantile", "IDB")
    return broadcast(
        (delta, kappa, gamma_param, u) -> _continuous_quantile_scalar(t -> _idb_pdf(delta, kappa, gamma_param, t), 0.0, Inf, u),
        d.delta,
        d.kappa,
        d.gamma_param,
        u,
    )
end

rand(rng::AbstractRNG, d::IDB) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::IDB, n::Integer) = quantile(d, rand(rng, n))
