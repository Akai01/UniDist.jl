struct Laplace{Ta,Tb} <: ContinuousDistribution
    alpha1::Ta
    alpha2::Tb
end

support(::Laplace) = (-Inf, Inf)

function _laplace_pdf(alpha1::Real, alpha2::Real, x::Real)
    if x < 0
        return exp(x / alpha1) / (alpha1 + alpha2)
    end
    return exp(-x / alpha2) / (alpha1 + alpha2)
end

pdf(d::Laplace, x) = broadcast((a1, a2, x) -> _laplace_pdf(a1, a2, x), d.alpha1, d.alpha2, x)

function cdf(d::Laplace, x)
    return broadcast((a1, a2, x) -> begin
        if x < 0
            (a1 / (a1 + a2)) * exp(x / a1)
        else
            1 - (a2 / (a1 + a2)) * exp(-x / a2)
        end
    end, d.alpha1, d.alpha2, x)
end

function quantile(d::Laplace, u)
    return broadcast((a1, a2, u) -> begin
        if u <= 0
            -Inf
        elseif u >= 1
            Inf
        else
            p = a1 / (a1 + a2)
            if u < p
                a1 * log(u / p)
            else
                -a2 * log((1 - u) / (1 - p))
            end
        end
    end, d.alpha1, d.alpha2, u)
end

rand(rng::AbstractRNG, d::Laplace) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::Laplace, n::Integer) = quantile(d, rand(rng, n))

mean(d::Laplace) = @. d.alpha2 - d.alpha1
var(d::Laplace) = @. d.alpha1 ^ 2 + d.alpha2 ^ 2
