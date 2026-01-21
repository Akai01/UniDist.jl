struct LogisticExponential{Ta,Tb} <: ContinuousDistribution
    alpha::Ta
    beta::Tb
end

support(::LogisticExponential) = (0.0, Inf)

function _logistic_exponential_pdf(alpha::Real, beta::Real, x::Real)
    if x <= 0
        return 0.0
    end
    ex = exp(alpha * x) - 1
    num = alpha * beta * ex ^ (beta - 1) * exp(alpha * x)
    den = (1 + ex ^ beta) ^ 2
    return num / den
end

pdf(d::LogisticExponential, x) = broadcast((alpha, beta, x) -> _logistic_exponential_pdf(alpha, beta, x), d.alpha, d.beta, x)

cdf(d::LogisticExponential, x) = broadcast((alpha, beta, x) -> x <= 0 ? 0.0 : (exp(alpha * x) - 1) ^ beta / (1 + (exp(alpha * x) - 1) ^ beta), d.alpha, d.beta, x)

quantile(d::LogisticExponential, u) = broadcast((alpha, beta, u) -> begin
    if u <= 0
        0.0
    elseif u >= 1
        Inf
    else
        log(1 + (u / (1 - u)) ^ (1 / beta)) / alpha
    end
end, d.alpha, d.beta, u)

rand(rng::AbstractRNG, d::LogisticExponential) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::LogisticExponential, n::Integer) = quantile(d, rand(rng, n))
