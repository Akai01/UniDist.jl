import SpecialFunctions: gamma

struct Erlang{Talpha,Tn} <: ContinuousDistribution
    alpha::Talpha
    n::Tn
end

support(::Erlang) = (0.0, Inf)

function _erlang_pdf(alpha::Real, n::Real, x::Real)
    if x <= 0
        return 0.0
    end
    return x ^ (n - 1) * exp(-x / alpha) / (alpha ^ n * gamma(n))
end

pdf(d::Erlang, x) = broadcast((alpha, n, x) -> _erlang_pdf(alpha, n, x), d.alpha, d.n, x)

function cdf(d::Erlang, x)
    return broadcast(
        (alpha, n, x) -> _continuous_cdf_scalar(t -> _erlang_pdf(alpha, n, t), 0.0, Inf, x),
        d.alpha,
        d.n,
        x,
    )
end

function quantile(d::Erlang, u)
    _info_numeric("quantile", "Erlang")
    return broadcast(
        (alpha, n, u) -> _continuous_quantile_scalar(t -> _erlang_pdf(alpha, n, t), 0.0, Inf, u),
        d.alpha,
        d.n,
        u,
    )
end

rand(rng::AbstractRNG, d::Erlang) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::Erlang, n::Integer) = quantile(d, rand(rng, n))

mean(d::Erlang) = @. d.alpha * d.n
var(d::Erlang) = @. d.alpha ^ 2 * d.n
