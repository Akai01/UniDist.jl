struct Rectangular{T} <: DiscreteDistribution
    n::T
end

function support(d::Rectangular)
    return d.n isa AbstractArray ? nothing : (0:Int(d.n))
end

function _rectangular_pmf(n::Real, x::Real)
    if !is_integer_value(x) || x < 0 || x > n
        return 0.0
    end
    return 1 / (n + 1)
end

pdf(d::Rectangular, x) = broadcast((n, x) -> _rectangular_pmf(n, x), d.n, x)

function cdf(d::Rectangular, x)
    return broadcast(
        (n, x) -> x < 0 ? 0.0 : (x >= n ? 1.0 : (floor(Int, x) + 1) / (n + 1)),
        d.n,
        x,
    )
end

function quantile(d::Rectangular, u)
    return broadcast(
        (n, u) -> begin
            if u <= 0
                return 0
            end
            if u >= 1
                return Int(n)
            end
            return floor(Int, u * (n + 1))
        end,
        d.n,
        u,
    )
end

function rand(rng::AbstractRNG, d::Rectangular)
    shape = broadcast_shape(d.n)
    if isempty(shape)
        return quantile(d, rand(rng))
    end
    u = rand(rng, shape)
    return quantile(d, u)
end

function rand(rng::AbstractRNG, d::Rectangular, n::Integer)
    shape = broadcast_shape(d.n)
    if isempty(shape)
        u = rand(rng, n)
        return quantile(d, u)
    end
    u = rand(rng, (n, shape...))
    return quantile(d, u)
end

mean(d::Rectangular) = @. d.n / 2
var(d::Rectangular) = @. d.n * (d.n + 2) / 12
