"""
    hdi(dist, mass=0.9; grid=2048, eps=1e-6)

A practical HDI/HPD approximation.
- Continuous: grid-based density scan on [q(eps), q(1-eps)].
- Discrete: exact-ish via sorting pmf over finite support.

Currently expects scalar distribution parameters.
"""
function hdi end

function hdi(dist::ContinuousDistribution, mass::Real=0.9; grid::Int=2048, eps::Real=1e-6)
    a = quantile(dist, eps)
    b = quantile(dist, 1 - eps)
    if a isa AbstractArray || b isa AbstractArray
        throw(ArgumentError("hdi currently supports scalar parameters"))
    end
    xs = range(a, b; length=grid)
    fs = pdf(dist, xs)
    ord = sortperm(fs; rev=true)
    w = (b - a) / (grid - 1)
    cum = 0.0
    keep = falses(grid)
    for i in ord
        keep[i] = true
        cum += fs[i] * w
        cum >= mass && break
    end
    idx = findall(keep)
    return (minimum(xs[idx]), maximum(xs[idx]))
end

function hdi(dist::DiscreteDistribution, mass::Real=0.9)
    supp = support(dist)
    supp === nothing && throw(ArgumentError("support required for discrete hdi"))
    xs = collect(supp)
    fs = pdf(dist, xs)
    if fs isa AbstractArray && length(size(fs)) > 1
        throw(ArgumentError("hdi currently supports scalar parameters"))
    end
    ord = sortperm(fs; rev=true)
    cum = 0.0
    keep = falses(length(xs))
    for i in ord
        keep[i] = true
        cum += fs[i]
        cum >= mass && break
    end
    idx = findall(keep)
    return (minimum(xs[idx]), maximum(xs[idx]))
end
