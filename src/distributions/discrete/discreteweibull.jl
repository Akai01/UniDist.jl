struct DiscreteWeibull{Tp,Tb} <: DiscreteDistribution
    p::Tp
    beta::Tb
end

support(::DiscreteWeibull) = nothing

function _discrete_weibull_pmf(p::Real, beta::Real, x::Real)
    if !is_integer_value(x) || x < 0
        return 0.0
    end
    k = Int(x)
    logq = log1p(-p)
    a = exp(logq * (k ^ beta))
    b = exp(logq * ((k + 1) ^ beta))
    return a - b
end

pdf(d::DiscreteWeibull, x) =
    broadcast((p, beta, x) -> _discrete_weibull_pmf(p, beta, x), d.p, d.beta, x)

function cdf(d::DiscreteWeibull, x)
    return broadcast(
        (p, beta, x) -> begin
            if x < 0
                0.0
            else
                1 - exp(log1p(-p) * ((floor(Int, x) + 1) ^ beta))
            end
        end,
        d.p,
        d.beta,
        x,
    )
end

function quantile(d::DiscreteWeibull, u)
    return broadcast(
        (p, beta, u) -> _discrete_quantile_scalar(
            k -> _discrete_weibull_pmf(p, beta, k),
            0,
            nothing,
            u,
        ),
        d.p,
        d.beta,
        u,
    )
end

function rand(rng::AbstractRNG, d::DiscreteWeibull)
    shape = broadcast_shape(d.p, d.beta)
    if isempty(shape)
        return quantile(d, rand(rng))
    end
    u = rand(rng, shape)
    return quantile(d, u)
end

function rand(rng::AbstractRNG, d::DiscreteWeibull, n::Integer)
    shape = broadcast_shape(d.p, d.beta)
    if isempty(shape)
        u = rand(rng, n)
        return quantile(d, u)
    end
    u = rand(rng, (n, shape...))
    return quantile(d, u)
end
