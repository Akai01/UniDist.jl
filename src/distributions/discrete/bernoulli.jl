struct Bernoulli{T} <: DiscreteDistribution
    p::T
end

support(::Bernoulli) = 0:1

function _bernoulli_pmf(p::Real, x::Real)
    if x == 0
        return 1 - p
    elseif x == 1
        return p
    end
    return 0.0
end

pdf(d::Bernoulli, x) = broadcast((p, x) -> _bernoulli_pmf(p, x), d.p, x)

function cdf(d::Bernoulli, x)
    return broadcast((p, x) -> x < 0 ? 0.0 : (x < 1 ? 1 - p : 1.0), d.p, x)
end

quantile(d::Bernoulli, u) = broadcast((p, u) -> u < 1 - p ? 0 : 1, d.p, u)

function rand(rng::AbstractRNG, d::Bernoulli)
    shape = broadcast_shape(d.p)
    if isempty(shape)
        return rand(rng) < d.p ? 1 : 0
    end
    u = rand(rng, shape)
    return ifelse.(u .< d.p, 1, 0)
end

function rand(rng::AbstractRNG, d::Bernoulli, n::Integer)
    shape = broadcast_shape(d.p)
    if isempty(shape)
        u = rand(rng, n)
        return ifelse.(u .< d.p, 1, 0)
    end
    u = rand(rng, (n, shape...))
    return ifelse.(u .< d.p, 1, 0)
end

mean(d::Bernoulli) = d.p
var(d::Bernoulli) = @. d.p * (1 - d.p)
