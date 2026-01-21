struct StandardTriangular <: ContinuousDistribution end

support(::StandardTriangular) = (-1.0, 1.0)

function _standard_triangular_pdf(x::Real)
    if x <= -1 || x >= 1
        return 0.0
    elseif x < 0
        return x + 1
    else
        return 1 - x
    end
end

pdf(::StandardTriangular, x) = @. _standard_triangular_pdf(x)

function cdf(::StandardTriangular, x)
    return @. ifelse(x <= -1, 0.0, ifelse(x < 0, 0.5 * x ^ 2 + x + 0.5, ifelse(x < 1, -0.5 * x ^ 2 + x + 0.5, 1.0)))
end

function quantile(::StandardTriangular, u)
    return @. ifelse(u <= 0, -1.0, ifelse(u >= 1, 1.0, ifelse(u < 0.5, -1 + sqrt(2 * u), 1 - sqrt(2 * (1 - u)))))
end

rand(rng::AbstractRNG, d::StandardTriangular) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::StandardTriangular, n::Integer) = quantile(d, rand(rng, n))

mean(::StandardTriangular) = 0.0
var(::StandardTriangular) = 1 / 6
