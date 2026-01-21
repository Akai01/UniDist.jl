struct Hypoexponential{Ta} <: ContinuousDistribution
    alpha::Ta
end

support(::Hypoexponential) = (0.0, Inf)

function _hypoexponential_pdf(alpha::AbstractVector, x::Real)
    if x <= 0
        return 0.0
    end
    s = 0.0
    n = length(alpha)
    for i in 1:n
        term = exp(-x / alpha[i]) / alpha[i]
        for j in 1:n
            if j != i
                term *= alpha[i] / (alpha[i] - alpha[j])
            end
        end
        s += term
    end
    return s
end

pdf(d::Hypoexponential, x) = broadcast(x -> _hypoexponential_pdf(d.alpha, x), x)

function cdf(d::Hypoexponential, x)
    _warn_numeric("cdf", "Hypoexponential")
    return broadcast(x -> _continuous_cdf_scalar(t -> _hypoexponential_pdf(d.alpha, t), 0.0, Inf, x), x)
end

function quantile(d::Hypoexponential, u)
    _warn_numeric("quantile", "Hypoexponential")
    return broadcast(u -> _continuous_quantile_scalar(t -> _hypoexponential_pdf(d.alpha, t), 0.0, Inf, u), u)
end

rand(rng::AbstractRNG, d::Hypoexponential) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::Hypoexponential, n::Integer) = quantile(d, rand(rng, n))
