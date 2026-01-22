import SpecialFunctions: gamma

struct Chi{T} <: ContinuousDistribution
    n::T
end

support(::Chi) = (0.0, Inf)

function _chi_pdf(n::Real, x::Real)
    if x <= 0
        return 0.0
    end
    return (1 / (2 ^ (n / 2 - 1) * gamma(n / 2))) * x ^ (n - 1) * exp(-x ^ 2 / 2)
end

pdf(d::Chi, x) = broadcast((n, x) -> _chi_pdf(n, x), d.n, x)

function cdf(d::Chi, x)
    _info_numeric("cdf", "Chi")
    return broadcast((n, x) -> _continuous_cdf_scalar(t -> _chi_pdf(n, t), 0.0, Inf, x), d.n, x)
end

function quantile(d::Chi, u)
    _info_numeric("quantile", "Chi")
    return broadcast((n, u) -> _continuous_quantile_scalar(t -> _chi_pdf(n, t), 0.0, Inf, u), d.n, u)
end

rand(rng::AbstractRNG, d::Chi) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::Chi, n::Integer) = quantile(d, rand(rng, n))

mean(d::Chi) = @. sqrt(2) * gamma((d.n + 1) / 2) / gamma(d.n / 2)
var(d::Chi) = @. d.n - (sqrt(2) * gamma((d.n + 1) / 2) / gamma(d.n / 2)) ^ 2
