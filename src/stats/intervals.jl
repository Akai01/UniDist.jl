"""
    interval(dist, alpha=0.05)

Equal-tailed interval [q(alpha/2), q(1-alpha/2)].
"""
function interval(dist::AbstractDistribution, alpha::Real=0.05)
    lo = quantile(dist, alpha / 2)
    hi = quantile(dist, 1 - alpha / 2)
    return (lo, hi)
end
