import SpecialFunctions: zeta

struct Zeta{T} <: DiscreteDistribution
    alpha::T
end

support(::Zeta) = nothing

function _zeta_pmf(alpha::Real, x::Real)
    if !is_integer_value(x) || x < 1
        return 0.0
    end
    k = Int(x)
    return 1 / (k ^ alpha * zeta(alpha))
end

pdf(d::Zeta, x) = broadcast((alpha, x) -> _zeta_pmf(alpha, x), d.alpha, x)

function cdf(d::Zeta, x)
    return broadcast(
        (alpha, x) -> _discrete_cdf_scalar(k -> _zeta_pmf(alpha, k), 1, nothing, x),
        d.alpha,
        x,
    )
end

function quantile(d::Zeta, u)
    return broadcast(
        (alpha, u) -> _discrete_quantile_scalar(k -> _zeta_pmf(alpha, k), 1, nothing, u),
        d.alpha,
        u,
    )
end

function rand(rng::AbstractRNG, d::Zeta)
    shape = broadcast_shape(d.alpha)
    if isempty(shape)
        return quantile(d, rand(rng))
    end
    u = rand(rng, shape)
    return quantile(d, u)
end

function rand(rng::AbstractRNG, d::Zeta, n::Integer)
    shape = broadcast_shape(d.alpha)
    if isempty(shape)
        u = rand(rng, n)
        return quantile(d, u)
    end
    u = rand(rng, (n, shape...))
    return quantile(d, u)
end

mean(d::Zeta) = @. zeta(d.alpha - 1) / zeta(d.alpha)
var(d::Zeta) = @. zeta(d.alpha - 2) / zeta(d.alpha) - (zeta(d.alpha - 1) / zeta(d.alpha)) ^ 2
