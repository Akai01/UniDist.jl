struct Minimax{Tbeta,Tgamma} <: ContinuousDistribution
    beta::Tbeta
    gamma::Tgamma
end

support(::Minimax) = (0.0, 1.0)

function _minimax_pdf(beta::Real, gamma::Real, x::Real)
    if x <= 0 || x >= 1
        return 0.0
    end
    return beta * gamma * x ^ (beta - 1) * (1 - x ^ beta) ^ (gamma - 1)
end

pdf(d::Minimax, x) = broadcast((beta, gamma, x) -> _minimax_pdf(beta, gamma, x), d.beta, d.gamma, x)

cdf(d::Minimax, x) = broadcast((beta, gamma, x) -> x <= 0 ? 0.0 : (x >= 1 ? 1.0 : 1 - (1 - x ^ beta) ^ gamma), d.beta, d.gamma, x)

quantile(d::Minimax, u) = broadcast((beta, gamma, u) -> begin
    if u <= 0
        0.0
    elseif u >= 1
        1.0
    else
        (1 - (1 - u) ^ (1 / gamma)) ^ (1 / beta)
    end
end, d.beta, d.gamma, u)

rand(rng::AbstractRNG, d::Minimax) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::Minimax, n::Integer) = quantile(d, rand(rng, n))
