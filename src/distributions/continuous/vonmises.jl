import SpecialFunctions: besseli

struct VonMises{Tkappa,Tmu} <: ContinuousDistribution
    kappa::Tkappa
    mu::Tmu
end

support(::VonMises) = (0.0, 2 * pi)

function _vonmises_pdf(kappa::Real, mu::Real, x::Real)
    if x <= 0 || x >= 2 * pi
        return 0.0
    end
    return exp(kappa * cos(x - mu)) / (2 * pi * besseli(0, kappa))
end

pdf(d::VonMises, x) = broadcast((kappa, mu, x) -> _vonmises_pdf(kappa, mu, x), d.kappa, d.mu, x)

function cdf(d::VonMises, x)
    _warn_numeric("cdf", "VonMises")
    return broadcast(
        (kappa, mu, x) -> _continuous_cdf_scalar(t -> _vonmises_pdf(kappa, mu, t), 0.0, 2 * pi, x),
        d.kappa,
        d.mu,
        x,
    )
end

function quantile(d::VonMises, u)
    _warn_numeric("quantile", "VonMises")
    return broadcast(
        (kappa, mu, u) -> _continuous_quantile_scalar(t -> _vonmises_pdf(kappa, mu, t), 0.0, 2 * pi, u),
        d.kappa,
        d.mu,
        u,
    )
end

rand(rng::AbstractRNG, d::VonMises) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::VonMises, n::Integer) = quantile(d, rand(rng, n))
