import QuadGK: quadgk
import Logging: @warn

const _SERIES_TERMS = 100

function _continuous_cdf_scalar(pdf::Function, a::Real, b::Real, x::Real)
    if x <= a
        return 0.0
    end
    if x >= b
        return 1.0
    end
    val, _ = quadgk(pdf, a, x)
    return clamp(val, 0.0, 1.0)
end

function _bracket_quantile(cdf::Function, a::Real, b::Real, u::Real; max_expand::Int=200)
    lo = isfinite(a) ? a : -1.0
    hi = isfinite(b) ? b : 1.0
    if lo >= hi
        hi = lo + 1.0
    end

    if !isfinite(a)
        step = 1.0
        for _ in 1:max_expand
            if cdf(lo) <= u
                break
            end
            hi = lo
            lo -= step
            step *= 2.0
        end
    end

    if !isfinite(b)
        step = 1.0
        for _ in 1:max_expand
            if cdf(hi) >= u
                break
            end
            lo = hi
            hi += step
            step *= 2.0
        end
    end

    return lo, hi
end

function _continuous_quantile_scalar(
    pdf::Function,
    a::Real,
    b::Real,
    u::Real;
    tol::Real=1e-8,
    max_iter::Int=200,
)
    if u <= 0
        return a
    end
    if u >= 1
        return b
    end
    cdf = x -> _continuous_cdf_scalar(pdf, a, b, x)
    lo, hi = _bracket_quantile(cdf, a, b, u)
    for _ in 1:max_iter
        mid = (lo + hi) / 2
        cmid = cdf(mid)
        if abs(cmid - u) <= tol
            return mid
        end
        if cmid < u
            lo = mid
        else
            hi = mid
        end
    end
    return (lo + hi) / 2
end

function _info_numeric(method::AbstractString, distname::AbstractString)
    @info("$method for $distname uses numerical integration/bisection")
end
