const FrameFunWrapper = FunctionWrappersWrapper

@inline function _build_frame_fun_wrapper(fun::Function, argtypes, rettypes)
    return FunctionWrappersWrapper(fun, argtypes, rettypes; cache=NoCache(), policy=AllowAll())
end

@inline function _frame_point_fun_wrapper(::Val{O}, ::Type{T}, fun::Function) where {O,T}
    return _build_frame_fun_wrapper(fun, (Tuple{T},), (Translation{O,T},))
end

@inline function _frame_axes_fun_wrapper(::Val{O}, ::Type{T}, fun::Function) where {O,T}
    return _build_frame_fun_wrapper(fun, (Tuple{T},), (Rotation{O,T},))
end
