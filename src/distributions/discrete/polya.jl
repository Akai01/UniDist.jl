import SpecialFunctions: loggamma

struct Polya{Tn,Tp,Tb} <: DiscreteDistribution
    n::Tn
    p::Tp
    beta::Tb
end

function support(d::Polya)
    return d.n isa AbstractArray ? nothing : (0:Int(d.n))
end

function _polya_pmf(n::Real, p::Real, beta::Real, x::Real)
    if !is_integer_value(x) || x < 0 || x > n
        return 0.0
    end
    a = p / beta
    b = (1 - p) / beta
    k = Int(x)
    logpmf = loggamma(n + 1) - loggamma(k + 1) - loggamma(n - k + 1) +
        loggamma(k + a) + loggamma(n - k + b) + loggamma(a + b) -
        loggamma(n + a + b) - loggamma(a) - loggamma(b)
    return exp(logpmf)
end

pdf(d::Polya, x) = broadcast((n, p, beta, x) -> _polya_pmf(n, p, beta, x), d.n, d.p, d.beta, x)

function cdf(d::Polya, x)
    return broadcast(
        (n, p, beta, x) -> _discrete_cdf_scalar(k -> _polya_pmf(n, p, beta, k), 0, Int(n), x),
        d.n,
        d.p,
        d.beta,
        x,
    )
end

function quantile(d::Polya, u)
    return broadcast(
        (n, p, beta, u) -> _discrete_quantile_scalar(k -> _polya_pmf(n, p, beta, k), 0, Int(n), u),
        d.n,
        d.p,
        d.beta,
        u,
    )
end

function rand(rng::AbstractRNG, d::Polya)
    shape = broadcast_shape(d.n, d.p, d.beta)
    if isempty(shape)
        return quantile(d, rand(rng))
    end
    u = rand(rng, shape)
    return quantile(d, u)
end

function rand(rng::AbstractRNG, d::Polya, n::Integer)
    shape = broadcast_shape(d.n, d.p, d.beta)
    if isempty(shape)
        u = rand(rng, n)
        return quantile(d, u)
    end
    u = rand(rng, (n, shape...))
    return quantile(d, u)
end

mean(d::Polya) = @. d.n * d.p

function var(d::Polya)
    a = @. d.p / d.beta
    b = @. (1 - d.p) / d.beta
    return @. d.n * a * b * (a + b + d.n) / ((a + b) ^ 2 * (a + b + 1))
end
