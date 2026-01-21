struct Logarithm{T} <: DiscreteDistribution
    c::T
end

support(::Logarithm) = nothing

function _logarithm_pmf(c::Real, x::Real)
    if !is_integer_value(x) || x < 1
        return 0.0
    end
    k = Int(x)
    return -((1 - c) ^ k) / (k * log(c))
end

pdf(d::Logarithm, x) = broadcast((c, x) -> _logarithm_pmf(c, x), d.c, x)

function cdf(d::Logarithm, x)
    return broadcast(
        (c, x) -> _discrete_cdf_scalar(k -> _logarithm_pmf(c, k), 1, nothing, x),
        d.c,
        x,
    )
end

function quantile(d::Logarithm, u)
    return broadcast(
        (c, u) -> _discrete_quantile_scalar(k -> _logarithm_pmf(c, k), 1, nothing, u),
        d.c,
        u,
    )
end

function rand(rng::AbstractRNG, d::Logarithm)
    shape = broadcast_shape(d.c)
    if isempty(shape)
        return quantile(d, rand(rng))
    end
    u = rand(rng, shape)
    return quantile(d, u)
end

function rand(rng::AbstractRNG, d::Logarithm, n::Integer)
    shape = broadcast_shape(d.c)
    if isempty(shape)
        u = rand(rng, n)
        return quantile(d, u)
    end
    u = rand(rng, (n, shape...))
    return quantile(d, u)
end

mean(d::Logarithm) = @. (d.c - 1) / (d.c * log(d.c))

function var(d::Logarithm)
    return @. (d.c - 1) * (log(d.c) + 1 - d.c) / ((log(d.c)) ^ 2 * (d.c ^ 2))
end
