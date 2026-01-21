import SpecialFunctions: loggamma

struct GammaPoisson{Ta,Tb} <: DiscreteDistribution
    alpha::Ta
    beta::Tb
end

support(::GammaPoisson) = nothing

function _gamma_poisson_pmf(alpha::Real, beta::Real, x::Real)
    if !is_integer_value(x) || x < 0
        return 0.0
    end
    k = Int(x)
    logpmf = loggamma(k + beta) + k * log(alpha) - loggamma(beta) -
        (beta + k) * log(1 + alpha) - loggamma(k + 1)
    return exp(logpmf)
end

pdf(d::GammaPoisson, x) =
    broadcast((alpha, beta, x) -> _gamma_poisson_pmf(alpha, beta, x), d.alpha, d.beta, x)

function cdf(d::GammaPoisson, x)
    return broadcast(
        (alpha, beta, x) -> _discrete_cdf_scalar(k -> _gamma_poisson_pmf(alpha, beta, k), 0, nothing, x),
        d.alpha,
        d.beta,
        x,
    )
end

function quantile(d::GammaPoisson, u)
    return broadcast(
        (alpha, beta, u) -> _discrete_quantile_scalar(k -> _gamma_poisson_pmf(alpha, beta, k), 0, nothing, u),
        d.alpha,
        d.beta,
        u,
    )
end

function rand(rng::AbstractRNG, d::GammaPoisson)
    shape = broadcast_shape(d.alpha, d.beta)
    if isempty(shape)
        return quantile(d, rand(rng))
    end
    u = rand(rng, shape)
    return quantile(d, u)
end

function rand(rng::AbstractRNG, d::GammaPoisson, n::Integer)
    shape = broadcast_shape(d.alpha, d.beta)
    if isempty(shape)
        u = rand(rng, n)
        return quantile(d, u)
    end
    u = rand(rng, (n, shape...))
    return quantile(d, u)
end

mean(d::GammaPoisson) = @. d.alpha * d.beta
var(d::GammaPoisson) = @. d.alpha * d.beta * (1 + d.alpha)
