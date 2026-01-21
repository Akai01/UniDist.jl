struct Binomial{Tn,Tp} <: DiscreteDistribution
    n::Tn
    p::Tp
end

function support(d::Binomial)
    return d.n isa AbstractArray ? nothing : (0:Int(d.n))
end

function _binomial_pmf(n::Real, p::Real, x::Real)
    if !is_integer_value(x) || x < 0 || x > n
        return 0.0
    end
    k = Int(x)
    return binomial(Int(n), k) * p^k * (1 - p)^(Int(n) - k)
end

pdf(d::Binomial, x) = broadcast((n, p, x) -> _binomial_pmf(n, p, x), d.n, d.p, x)

function cdf(d::Binomial, x)
    return broadcast(
        (n, p, x) -> _discrete_cdf_scalar(k -> _binomial_pmf(n, p, k), 0, Int(n), x),
        d.n,
        d.p,
        x,
    )
end

function quantile(d::Binomial, u)
    return broadcast(
        (n, p, u) -> _discrete_quantile_scalar(k -> _binomial_pmf(n, p, k), 0, Int(n), u),
        d.n,
        d.p,
        u,
    )
end

function rand(rng::AbstractRNG, d::Binomial)
    shape = broadcast_shape(d.n, d.p)
    if isempty(shape)
        return quantile(d, rand(rng))
    end
    u = rand(rng, shape)
    return quantile(d, u)
end

function rand(rng::AbstractRNG, d::Binomial, n::Integer)
    shape = broadcast_shape(d.n, d.p)
    if isempty(shape)
        u = rand(rng, n)
        return quantile(d, u)
    end
    u = rand(rng, (n, shape...))
    return quantile(d, u)
end

mean(d::Binomial) = @. d.n * d.p
var(d::Binomial) = @. d.n * d.p * (1 - d.p)
