import SpecialFunctions: gamma

struct F{Tn1,Tn2} <: ContinuousDistribution
    n1::Tn1
    n2::Tn2
end

support(::F) = (0.0, Inf)

function _f_pdf(n1::Real, n2::Real, x::Real)
    if x <= 0
        return 0.0
    end
    num = gamma((n1 + n2) / 2) * (n1 / n2) ^ (n1 / 2) * x ^ (n1 / 2 - 1)
    den = gamma(n1 / 2) * gamma(n2 / 2) * (1 + (n1 / n2) * x) ^ ((n1 + n2) / 2)
    return num / den
end

pdf(d::F, x) = broadcast((n1, n2, x) -> _f_pdf(n1, n2, x), d.n1, d.n2, x)

function cdf(d::F, x)
    return broadcast((n1, n2, x) -> _continuous_cdf_scalar(t -> _f_pdf(n1, n2, t), 0.0, Inf, x), d.n1, d.n2, x)
end

function quantile(d::F, u)
    _info_numeric("quantile", "F")
    return broadcast((n1, n2, u) -> _continuous_quantile_scalar(t -> _f_pdf(n1, n2, t), 0.0, Inf, u), d.n1, d.n2, u)
end

rand(rng::AbstractRNG, d::F) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::F, n::Integer) = quantile(d, rand(rng, n))

mean(d::F) = @. d.n2 > 2 ? d.n2 / (d.n2 - 2) : Inf
var(d::F) = @. d.n2 > 4 ? 2 * d.n2 ^ 2 * (d.n1 + d.n2 - 2) / (d.n1 * (d.n2 - 2) ^ 2 * (d.n2 - 4)) : Inf
