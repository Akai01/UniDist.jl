struct Zipf{Ta,Tn} <: DiscreteDistribution
    alpha::Ta
    n::Tn
end

function support(d::Zipf)
    return d.n isa AbstractArray ? nothing : (1:Int(d.n))
end

function _zipf_norm(alpha::Real, n::Int)
    s = 0.0
    for i in 1:n
        s += i ^ (-alpha)
    end
    return s
end

function _zipf_pmf(alpha::Real, n::Real, x::Real)
    if !is_integer_value(x) || x < 1 || x > n
        return 0.0
    end
    k = Int(x)
    norm = _zipf_norm(alpha, Int(n))
    return 1 / (k ^ alpha * norm)
end

pdf(d::Zipf, x) = broadcast((alpha, n, x) -> _zipf_pmf(alpha, n, x), d.alpha, d.n, x)

function cdf(d::Zipf, x)
    return broadcast(
        (alpha, n, x) -> _discrete_cdf_scalar(k -> _zipf_pmf(alpha, n, k), 1, Int(n), x),
        d.alpha,
        d.n,
        x,
    )
end

function quantile(d::Zipf, u)
    return broadcast(
        (alpha, n, u) -> _discrete_quantile_scalar(k -> _zipf_pmf(alpha, n, k), 1, Int(n), u),
        d.alpha,
        d.n,
        u,
    )
end

function rand(rng::AbstractRNG, d::Zipf)
    shape = broadcast_shape(d.alpha, d.n)
    if isempty(shape)
        return quantile(d, rand(rng))
    end
    u = rand(rng, shape)
    return quantile(d, u)
end

function rand(rng::AbstractRNG, d::Zipf, n::Integer)
    shape = broadcast_shape(d.alpha, d.n)
    if isempty(shape)
        u = rand(rng, n)
        return quantile(d, u)
    end
    u = rand(rng, (n, shape...))
    return quantile(d, u)
end

function mean(d::Zipf)
    return broadcast((alpha, n) -> begin
        n_int = Int(n)
        norm = _zipf_norm(alpha, n_int)
        s = 0.0
        for i in 1:n_int
            s += i ^ (1 - alpha)
        end
        s / norm
    end, d.alpha, d.n)
end

function var(d::Zipf)
    return broadcast((alpha, n) -> begin
        n_int = Int(n)
        norm = _zipf_norm(alpha, n_int)
        s1 = 0.0
        s2 = 0.0
        for i in 1:n_int
            s1 += i ^ (1 - alpha)
            s2 += i ^ (2 - alpha)
        end
        mu = s1 / norm
        s2 / norm - mu ^ 2
    end, d.alpha, d.n)
end
