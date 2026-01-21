struct DiscreteUniform{Ta,Tb} <: DiscreteDistribution
    a::Ta
    b::Tb
end

function support(d::DiscreteUniform)
    if d.a isa AbstractArray || d.b isa AbstractArray
        return nothing
    end
    return Int(d.a):Int(d.b)
end

function _discrete_uniform_pmf(a::Real, b::Real, x::Real)
    if !is_integer_value(x) || x < a || x > b
        return 0.0
    end
    return 1 / (b - a + 1)
end

pdf(d::DiscreteUniform, x) = broadcast((a, b, x) -> _discrete_uniform_pmf(a, b, x), d.a, d.b, x)

function cdf(d::DiscreteUniform, x)
    return broadcast(
        (a, b, x) -> x < a ? 0.0 : (x >= b ? 1.0 : (floor(Int, x) - a + 1) / (b - a + 1)),
        d.a,
        d.b,
        x,
    )
end

function quantile(d::DiscreteUniform, u)
    return broadcast(
        (a, b, u) -> begin
            if u <= 0
                return Int(a)
            end
            if u >= 1
                return Int(b)
            end
            return Int(a) + floor(Int, u * (b - a + 1))
        end,
        d.a,
        d.b,
        u,
    )
end

function rand(rng::AbstractRNG, d::DiscreteUniform)
    shape = broadcast_shape(d.a, d.b)
    if isempty(shape)
        return quantile(d, rand(rng))
    end
    u = rand(rng, shape)
    return quantile(d, u)
end

function rand(rng::AbstractRNG, d::DiscreteUniform, n::Integer)
    shape = broadcast_shape(d.a, d.b)
    if isempty(shape)
        u = rand(rng, n)
        return quantile(d, u)
    end
    u = rand(rng, (n, shape...))
    return quantile(d, u)
end

mean(d::DiscreteUniform) = @. (d.a + d.b) / 2
var(d::DiscreteUniform) = @. ((d.b - d.a + 1) ^ 2 - 1) / 12
