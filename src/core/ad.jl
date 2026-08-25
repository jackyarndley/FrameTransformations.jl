const FrameFunWrapper = FunctionWrappersWrapper
const _AUTODIFF_BACKEND = DI.AutoForwardDiff()

@inline derivative1(function_object, time) =
    DI.derivative(function_object, _AUTODIFF_BACKEND, time)

@inline derivative2(function_object, time) =
    DI.second_derivative(function_object, _AUTODIFF_BACKEND, time)

@inline derivative3(function_object, time) = DI.derivative(
    epoch -> DI.second_derivative(function_object, _AUTODIFF_BACKEND, epoch),
    _AUTODIFF_BACKEND,
    time,
)

@inline _value_derivatives(function_object, time, ::Val{1}) = (function_object(time),)

@inline _value_derivatives(function_object, time, ::Val{2}) =
    DI.value_and_derivative(function_object, _AUTODIFF_BACKEND, time)

@inline _value_derivatives(function_object, time, ::Val{3}) =
    DI.value_derivative_and_second_derivative(
        function_object, _AUTODIFF_BACKEND, time)

@inline function _value_derivatives(function_object, time, ::Val{4})
    lower = _value_derivatives(function_object, time, Val(3))
    return (lower..., derivative3(function_object, time))
end

@inline function _build_frame_fun_wrapper(fun, argtypes, rettypes)
    return FunctionWrappersWrapper(fun, argtypes, rettypes; cache=NoCache(), policy=AllowAll())
end

@inline function _frame_point_fun_wrapper(::Val{O}, ::Type{T}, fun) where {O,T}
    return _build_frame_fun_wrapper(fun, (Tuple{T},), (Translation{O,T},))
end

@inline function _frame_axes_fun_wrapper(::Val{O}, ::Type{T}, fun) where {O,T}
    return _build_frame_fun_wrapper(fun, (Tuple{T},), (Rotation{O,T},))
end

@inline function _frame_vector_fun_wrapper(::Val{O}, ::Type{T}, fun) where {O,T}
    return _build_frame_fun_wrapper(fun, (Tuple{T},), (SVector{3 * O,T},))
end
