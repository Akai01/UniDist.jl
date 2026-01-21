struct GammaNormal <: ContinuousDistribution end

support(::GammaNormal) = (-Inf, Inf)

function _gammanormal_error()
    throw(ArgumentError("Gamma-normal distribution is bivariate in the reference and is not implemented as a univariate distribution"))
end

pdf(::GammaNormal, x) = _gammanormal_error()
cdf(::GammaNormal, x) = _gammanormal_error()
quantile(::GammaNormal, u) = _gammanormal_error()

rand(rng::AbstractRNG, d::GammaNormal) = _gammanormal_error()
rand(rng::AbstractRNG, d::GammaNormal, n::Integer) = _gammanormal_error()
