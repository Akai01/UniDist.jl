struct Geometric{T} <: DiscreteDistribution
    p::T
end

support(::Geometric) = nothing

function _geometric_pmf(p::Real, x::Real)
    if !is_integer_value(x) || x < 0
        return 0.0
    end
    return p * (1 - p)^Int(x)
end

pdf(d::Geometric, x) = broadcast((p, x) -> _geometric_pmf(p, x), d.p, x)

function cdf(d::Geometric, x)
    return broadcast((p, x) -> x < 0 ? 0.0 : 1 - (1 - p)^(floor(Int, x) + 1), d.p, x)
end

function _geometric_quantile(p::Real, u::Real)
    if u <= 0
        return 0
    end
    if u >= 1
        return typemax(Int)
    end
    return max(0, ceil(Int, log1p(-u) / log1p(-p) - 1))
end

quantile(d::Geometric, u) = broadcast((p, u) -> _geometric_quantile(p, u), d.p, u)

function rand(rng::AbstractRNG, d::Geometric)
    shape = broadcast_shape(d.p)
    if isempty(shape)
        return _geometric_quantile(d.p, rand(rng))
    end
    u = rand(rng, shape)
    return quantile(d, u)
end

function rand(rng::AbstractRNG, d::Geometric, n::Integer)
    shape = broadcast_shape(d.p)
    if isempty(shape)
        u = rand(rng, n)
        return quantile(d, u)
    end
    u = rand(rng, (n, shape...))
    return quantile(d, u)
end

mean(d::Geometric) = @. (1 - d.p) / d.p
var(d::Geometric) = @. (1 - d.p) / (d.p ^ 2)
