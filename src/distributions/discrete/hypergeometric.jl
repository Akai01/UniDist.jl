struct Hypergeometric{T1,T2,T3} <: DiscreteDistribution
    n1::T1
    n2::T2
    n3::T3
end

function support(d::Hypergeometric)
    if d.n1 isa AbstractArray || d.n2 isa AbstractArray || d.n3 isa AbstractArray
        return nothing
    end
    lower = max(0, Int(d.n1 + d.n2 - d.n3))
    upper = min(Int(d.n1), Int(d.n2))
    return lower:upper
end

function _hypergeometric_pmf(n1::Real, n2::Real, n3::Real, x::Real)
    if !is_integer_value(x)
        return 0.0
    end
    k = Int(x)
    lower = max(0, Int(n1 + n2 - n3))
    upper = min(Int(n1), Int(n2))
    if k < lower || k > upper
        return 0.0
    end
    logpmf = _logchoose(n1, k) + _logchoose(n3 - n1, n2 - k) - _logchoose(n3, n2)
    return exp(logpmf)
end

pdf(d::Hypergeometric, x) =
    broadcast((n1, n2, n3, x) -> _hypergeometric_pmf(n1, n2, n3, x), d.n1, d.n2, d.n3, x)

function cdf(d::Hypergeometric, x)
    return broadcast(
        (n1, n2, n3, x) -> begin
            lower = max(0, Int(n1 + n2 - n3))
            upper = min(Int(n1), Int(n2))
            _discrete_cdf_scalar(k -> _hypergeometric_pmf(n1, n2, n3, k), lower, upper, x)
        end,
        d.n1,
        d.n2,
        d.n3,
        x,
    )
end

function quantile(d::Hypergeometric, u)
    return broadcast(
        (n1, n2, n3, u) -> begin
            lower = max(0, Int(n1 + n2 - n3))
            upper = min(Int(n1), Int(n2))
            _discrete_quantile_scalar(k -> _hypergeometric_pmf(n1, n2, n3, k), lower, upper, u)
        end,
        d.n1,
        d.n2,
        d.n3,
        u,
    )
end

function rand(rng::AbstractRNG, d::Hypergeometric)
    shape = broadcast_shape(d.n1, d.n2, d.n3)
    if isempty(shape)
        return quantile(d, rand(rng))
    end
    u = rand(rng, shape)
    return quantile(d, u)
end

function rand(rng::AbstractRNG, d::Hypergeometric, n::Integer)
    shape = broadcast_shape(d.n1, d.n2, d.n3)
    if isempty(shape)
        u = rand(rng, n)
        return quantile(d, u)
    end
    u = rand(rng, (n, shape...))
    return quantile(d, u)
end

mean(d::Hypergeometric) = @. d.n2 * d.n1 / d.n3

function var(d::Hypergeometric)
    return @. d.n2 * (d.n1 / d.n3) * (1 - d.n1 / d.n3) * (d.n3 - d.n2) / (d.n3 - 1)
end
