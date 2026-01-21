struct NegativeHypergeometric{T1,T2,T3} <: DiscreteDistribution
    n1::T1
    n2::T2
    n3::T3
end

function support(d::NegativeHypergeometric)
    return (d.n2 isa AbstractArray) ? nothing : (0:Int(d.n2))
end

function _negative_hypergeometric_pmf(n1::Real, n2::Real, n3::Real, x::Real)
    if !is_integer_value(x) || x < 0 || x > n2
        return 0.0
    end
    k = Int(x)
    logpmf = _logchoose(n1 + k - 1, k) +
        _logchoose(n3 - n1 + n2 - k - 1, n2 - k) -
        _logchoose(n3 + n2 - 1, n2)
    return exp(logpmf)
end

pdf(d::NegativeHypergeometric, x) =
    broadcast((n1, n2, n3, x) -> _negative_hypergeometric_pmf(n1, n2, n3, x), d.n1, d.n2, d.n3, x)

function cdf(d::NegativeHypergeometric, x)
    return broadcast(
        (n1, n2, n3, x) -> _discrete_cdf_scalar(
            k -> _negative_hypergeometric_pmf(n1, n2, n3, k),
            0,
            Int(n2),
            x,
        ),
        d.n1,
        d.n2,
        d.n3,
        x,
    )
end

function quantile(d::NegativeHypergeometric, u)
    return broadcast(
        (n1, n2, n3, u) -> _discrete_quantile_scalar(
            k -> _negative_hypergeometric_pmf(n1, n2, n3, k),
            0,
            Int(n2),
            u,
        ),
        d.n1,
        d.n2,
        d.n3,
        u,
    )
end

function rand(rng::AbstractRNG, d::NegativeHypergeometric)
    shape = broadcast_shape(d.n1, d.n2, d.n3)
    if isempty(shape)
        return quantile(d, rand(rng))
    end
    u = rand(rng, shape)
    return quantile(d, u)
end

function rand(rng::AbstractRNG, d::NegativeHypergeometric, n::Integer)
    shape = broadcast_shape(d.n1, d.n2, d.n3)
    if isempty(shape)
        u = rand(rng, n)
        return quantile(d, u)
    end
    u = rand(rng, (n, shape...))
    return quantile(d, u)
end

mean(d::NegativeHypergeometric) = @. d.n1 * d.n2 / d.n3

function var(d::NegativeHypergeometric)
    return @. d.n1 * d.n2 * (d.n3 + d.n2) / (d.n3 * (d.n3 + 1)) * (1 - d.n1 / d.n3)
end
