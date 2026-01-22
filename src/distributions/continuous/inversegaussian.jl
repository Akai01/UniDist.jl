struct InverseGaussian{Tlambda,Tmu} <: ContinuousDistribution
    lambda::Tlambda
    mu::Tmu
end

support(::InverseGaussian) = (0.0, Inf)

function _inverse_gaussian_pdf(lambda::Real, mu::Real, x::Real)
    if x <= 0
        return 0.0
    end
    return sqrt(lambda / (2 * pi * x ^ 3)) * exp(-lambda * (x - mu) ^ 2 / (2 * x * mu ^ 2))
end

pdf(d::InverseGaussian, x) = broadcast((lambda, mu, x) -> _inverse_gaussian_pdf(lambda, mu, x), d.lambda, d.mu, x)

function cdf(d::InverseGaussian, x)
    _info_numeric("cdf", "InverseGaussian")
    return broadcast(
        (lambda, mu, x) -> _continuous_cdf_scalar(t -> _inverse_gaussian_pdf(lambda, mu, t), 0.0, Inf, x),
        d.lambda,
        d.mu,
        x,
    )
end

function quantile(d::InverseGaussian, u)
    _info_numeric("quantile", "InverseGaussian")
    return broadcast(
        (lambda, mu, u) -> _continuous_quantile_scalar(t -> _inverse_gaussian_pdf(lambda, mu, t), 0.0, Inf, u),
        d.lambda,
        d.mu,
        u,
    )
end

rand(rng::AbstractRNG, d::InverseGaussian) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::InverseGaussian, n::Integer) = quantile(d, rand(rng, n))

mean(d::InverseGaussian) = d.mu
var(d::InverseGaussian) = @. d.mu ^ 3 / d.lambda
