struct Benford <: DiscreteDistribution end

support(::Benford) = 1:9

function _benford_pmf(x::Real)
    if !is_integer_value(x) || x < 1 || x > 9
        return 0.0
    end
    return log10(1 + 1 / Int(x))
end

pdf(::Benford, x) = broadcast(x -> _benford_pmf(x), x)

function cdf(::Benford, x)
    return broadcast(x -> begin
        if x < 1
            0.0
        elseif x >= 9
            1.0
        else
            log10(1 + floor(Int, x))
        end
    end, x)
end

function quantile(::Benford, u)
    return broadcast(u -> _discrete_quantile_scalar(_benford_pmf, 1, 9, u), u)
end

function rand(rng::AbstractRNG, d::Benford)
    return quantile(d, rand(rng))
end

function rand(rng::AbstractRNG, d::Benford, n::Integer)
    return quantile(d, rand(rng, n))
end

function mean(::Benford)
    s = 0.0
    for x in 1:9
        s += x * _benford_pmf(x)
    end
    return s
end

function var(::Benford)
    mu = mean(Benford())
    s = 0.0
    for x in 1:9
        s += (x ^ 2) * _benford_pmf(x)
    end
    return s - mu ^ 2
end
