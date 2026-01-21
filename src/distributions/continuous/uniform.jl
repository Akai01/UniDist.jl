struct Uniform{Ta,Tb} <: ContinuousDistribution
    a::Ta
    b::Tb
end

support(d::Uniform) = (d.a, d.b)

function _uniform_pdf(a::Real, b::Real, x::Real)
    if x <= a || x >= b
        return 0.0
    end
    return 1 / (b - a)
end

pdf(d::Uniform, x) = broadcast((a, b, x) -> _uniform_pdf(a, b, x), d.a, d.b, x)

function cdf(d::Uniform, x)
    return broadcast((a, b, x) -> x <= a ? 0.0 : (x >= b ? 1.0 : (x - a) / (b - a)), d.a, d.b, x)
end

function quantile(d::Uniform, u)
    return broadcast((a, b, u) -> u <= 0 ? a : (u >= 1 ? b : a + u * (b - a)), d.a, d.b, u)
end

rand(rng::AbstractRNG, d::Uniform) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::Uniform, n::Integer) = quantile(d, rand(rng, n))

mean(d::Uniform) = @. (d.a + d.b) / 2
var(d::Uniform) = @. (d.b - d.a) ^ 2 / 12
