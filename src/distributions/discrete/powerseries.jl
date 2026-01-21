struct PowerSeries{Ta,Tc} <: DiscreteDistribution
    a::Ta
    c::Tc
end

function support(d::PowerSeries)
    if d.a isa AbstractArray
        return 0:(length(d.a) - 1)
    end
    return nothing
end

function _powerseries_norm(a::AbstractVector, c::Real)
    s = 0.0
    for (idx, ax) in enumerate(a)
        x = idx - 1
        s += ax * c ^ x
    end
    return s
end

function _powerseries_pmf(a::AbstractVector, c::Real, x::Real)
    if !is_integer_value(x)
        return 0.0
    end
    k = Int(x)
    if k < 0 || k >= length(a)
        return 0.0
    end
    norm = _powerseries_norm(a, c)
    return a[k + 1] * c ^ k / norm
end

function pdf(d::PowerSeries, x)
    a = d.a
    return broadcast((c, x) -> _powerseries_pmf(a, c, x), d.c, x)
end

function cdf(d::PowerSeries, x)
    a = d.a
    return broadcast(
        (c, x) -> _discrete_cdf_scalar(k -> _powerseries_pmf(a, c, k), 0, length(a) - 1, x),
        d.c,
        x,
    )
end

function quantile(d::PowerSeries, u)
    a = d.a
    return broadcast(
        (c, u) -> _discrete_quantile_scalar(k -> _powerseries_pmf(a, c, k), 0, length(a) - 1, u),
        d.c,
        u,
    )
end

function rand(rng::AbstractRNG, d::PowerSeries)
    shape = broadcast_shape(d.c)
    if isempty(shape)
        return quantile(d, rand(rng))
    end
    u = rand(rng, shape)
    return quantile(d, u)
end

function rand(rng::AbstractRNG, d::PowerSeries, n::Integer)
    shape = broadcast_shape(d.c)
    if isempty(shape)
        u = rand(rng, n)
        return quantile(d, u)
    end
    u = rand(rng, (n, shape...))
    return quantile(d, u)
end

function mean(d::PowerSeries)
    a = d.a
    return broadcast((c) -> begin
        norm = _powerseries_norm(a, c)
        s = 0.0
        for (idx, ax) in enumerate(a)
            x = idx - 1
            s += x * ax * c ^ x
        end
        s / norm
    end, d.c)
end

function var(d::PowerSeries)
    a = d.a
    return broadcast((c) -> begin
        norm = _powerseries_norm(a, c)
        s1 = 0.0
        s2 = 0.0
        for (idx, ax) in enumerate(a)
            x = idx - 1
            w = ax * c ^ x
            s1 += x * w
            s2 += x ^ 2 * w
        end
        mu = s1 / norm
        s2 / norm - mu ^ 2
    end, d.c)
end
