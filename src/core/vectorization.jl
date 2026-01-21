param_shape(x::AbstractArray) = size(x)
param_shape(::Any) = ()

broadcast_shape(args...) = Base.broadcast_shape(map(param_shape, args)...)

is_integer_value(x::Real) = isfinite(x) && x == floor(x)
