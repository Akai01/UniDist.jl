struct Triangular{Ta,Tm,Tb} <: ContinuousDistribution
    a::Ta
    m::Tm
    b::Tb
end

support(d::Triangular) = (d.a, d.b)

function _triangular_pdf(a::Real, m::Real, b::Real, x::Real)
    if x <= a || x >= b
        return 0.0
    elseif x < m
        return 2 * (x - a) / ((b - a) * (m - a))
    else
        return 2 * (b - x) / ((b - a) * (b - m))
    end
end

pdf(d::Triangular, x) = broadcast((a, m, b, x) -> _triangular_pdf(a, m, b, x), d.a, d.m, d.b, x)

function cdf(d::Triangular, x)
    return broadcast((a, m, b, x) -> begin
        if x <= a
            0.0
        elseif x < m
            ((x - a) ^ 2) / ((b - a) * (m - a))
        elseif x < b
            1 - ((b - x) ^ 2) / ((b - a) * (b - m))
        else
            1.0
        end
    end, d.a, d.m, d.b, x)
end

function quantile(d::Triangular, u)
    return broadcast((a, m, b, u) -> begin
        if u <= 0
            a
        elseif u >= 1
            b
        else
            fm = (m - a) / (b - a)
            if u < fm
                a + sqrt(u * (b - a) * (m - a))
            else
                b - sqrt((1 - u) * (b - a) * (b - m))
            end
        end
    end, d.a, d.m, d.b, u)
end

rand(rng::AbstractRNG, d::Triangular) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::Triangular, n::Integer) = quantile(d, rand(rng, n))

mean(d::Triangular) = @. (d.a + d.m + d.b) / 3
var(d::Triangular) = @. (d.a ^ 2 + d.m ^ 2 + d.b ^ 2 - d.a * d.m - d.a * d.b - d.m * d.b) / 18
