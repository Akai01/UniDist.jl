struct Arctangent{Tlambda,Tphi} <: ContinuousDistribution
    lambda::Tlambda
    phi::Tphi
end

support(::Arctangent) = (0.0, Inf)

function _arctangent_pdf(lambda::Real, phi::Real, x::Real)
    if x < 0
        return 0.0
    end
    denom = (atan(lambda * phi) + 0.5 * pi) * (1 + (lambda ^ 2) * (x - phi) ^ 2)
    return lambda / denom
end

pdf(d::Arctangent, x) = broadcast((lambda, phi, x) -> _arctangent_pdf(lambda, phi, x), d.lambda, d.phi, x)

function cdf(d::Arctangent, x)
    return broadcast((lambda, phi, x) -> begin
        if x < 0
            0.0
        else
            num = atan(lambda * phi) - atan(-lambda * x + lambda * phi)
            den = atan(lambda * phi) * 2 + pi
            2 * num / den
        end
    end, d.lambda, d.phi, x)
end

function quantile(d::Arctangent, u)
    return broadcast((lambda, phi, u) -> begin
        if u <= 0
            0.0
        elseif u >= 1
            Inf
        else
            den = 2 * atan(lambda * phi) + pi
            theta = atan(lambda * phi) - 0.5 * u * den
            phi - tan(theta) / lambda
        end
    end, d.lambda, d.phi, u)
end

rand(rng::AbstractRNG, d::Arctangent) = quantile(d, rand(rng))
rand(rng::AbstractRNG, d::Arctangent, n::Integer) = quantile(d, rand(rng, n))
