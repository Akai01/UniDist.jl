struct KolmogorovSmirnov{T} <: ContinuousDistribution
    n::T
end

support(::KolmogorovSmirnov) = (0.0, 1.0)

function _ks_cdf_n1(x::Real)
    if x <= 0.5
        return 0.0
    elseif x < 1
        return 2 * x - 1
    else
        return 1.0
    end
end

function _ks_cdf_n2(x::Real)
    if x <= 0.25
        return 0.0
    elseif x < 0.5
        return (8 * x - 1) ^ 2 / 4
    elseif x < 1
        return 1 - 2 * (1 - x) ^ 2
    else
        return 1.0
    end
end

function _ks_cdf_asymptotic(n::Real, x::Real)
    if x <= 0
        return 0.0
    elseif x >= 1
        return 1.0
    end
    t = sqrt(n) * x
    s = 0.0
    for k in 1:_SERIES_TERMS
        s += (-1) ^ (k - 1) * exp(-2 * (k ^ 2) * (t ^ 2))
    end
    return clamp(1 - 2 * s, 0.0, 1.0)
end

function cdf(d::KolmogorovSmirnov, x)
    return broadcast((n, x) -> begin
        if x <= 0
            0.0
        elseif x >= 1
            1.0
        elseif n == 1
            _ks_cdf_n1(x)
        elseif n == 2
            _ks_cdf_n2(x)
        else
            _info_numeric("cdf", "KolmogorovSmirnov")
            _ks_cdf_asymptotic(n, x)
        end
    end, d.n, x)
end

function pdf(d::KolmogorovSmirnov, x)
    _info_numeric("pdf", "KolmogorovSmirnov")
    eps = 1e-6
    return broadcast((n, x) -> (cdf(KolmogorovSmirnov(n), x + eps) - cdf(KolmogorovSmirnov(n), x - eps)) / (2 * eps), d.n, x)
end

function quantile(d::KolmogorovSmirnov, u)
    _info_numeric("quantile", "KolmogorovSmirnov")
    return broadcast(
        (n, u) -> _continuous_quantile_scalar(t -> pdf(KolmogorovSmirnov(n), t), 0.0, 1.0, u),
        d.n,
        u,
    )
end

rand(rng::AbstractRNG, d::KolmogorovSmirnov) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::KolmogorovSmirnov, n::Integer) = quantile(d, rand(rng, n))
