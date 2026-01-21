import SpecialFunctions: loggamma

function _logchoose(n::Real, k::Real)
    if k < 0 || k > n
        return -Inf
    end
    return loggamma(n + 1) - loggamma(k + 1) - loggamma(n - k + 1)
end

function _discrete_cdf_scalar(pmf::Function, lower::Int, upper::Union{Int,Nothing}, x::Real)
    if !isfinite(x)
        return x < lower ? 0.0 : 1.0
    end
    if x < lower
        return 0.0
    end
    if upper !== nothing && x >= upper
        return 1.0
    end
    xmax = floor(Int, x)
    if upper !== nothing
        xmax = min(xmax, upper)
    end
    s = 0.0
    for k in lower:xmax
        s += pmf(k)
    end
    return s
end

function _discrete_quantile_scalar(
    pmf::Function,
    lower::Int,
    upper::Union{Int,Nothing},
    u::Real;
    max_iter::Int=10^7,
)
    if u <= 0
        return lower
    end
    if u >= 1
        return upper === nothing ? typemax(Int) : upper
    end
    s = 0.0
    k = lower
    while true
        s += pmf(k)
        if s >= u
            return k
        end
        if upper !== nothing && k >= upper
            return upper
        end
        k += 1
        if upper === nothing && (k - lower) >= max_iter
            return k
        end
    end
end
