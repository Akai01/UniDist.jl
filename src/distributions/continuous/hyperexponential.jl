struct Hyperexponential{Ta,Tp} <: ContinuousDistribution
    alpha::Ta
    p::Tp
end

support(::Hyperexponential) = (0.0, Inf)

function _hyperexponential_pdf(alpha::AbstractVector, p::AbstractVector, x::Real)
    if x <= 0
        return 0.0
    end
    s = 0.0
    for i in eachindex(alpha, p)
        s += p[i] * exp(-x / alpha[i]) / alpha[i]
    end
    return s
end

function _hyperexponential_cdf(alpha::AbstractVector, p::AbstractVector, x::Real)
    if x <= 0
        return 0.0
    end
    s = 0.0
    for i in eachindex(alpha, p)
        s += p[i] * exp(-x / alpha[i])
    end
    return 1 - s
end

pdf(d::Hyperexponential, x) = broadcast(x -> _hyperexponential_pdf(d.alpha, d.p, x), x)

cdf(d::Hyperexponential, x) = broadcast(x -> _hyperexponential_cdf(d.alpha, d.p, x), x)

function quantile(d::Hyperexponential, u)
    _info_numeric("quantile", "Hyperexponential")
    return broadcast(
        u -> _continuous_quantile_scalar(t -> _hyperexponential_pdf(d.alpha, d.p, t), 0.0, Inf, u),
        u,
    )
end

rand(rng::AbstractRNG, d::Hyperexponential) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::Hyperexponential, n::Integer) = quantile(d, rand(rng, n))
