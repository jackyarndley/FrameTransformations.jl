const TagAD1{T} = ForwardDiff.Tag{JSMDDiffTag,T}
const DualAD1{T} = ForwardDiff.Dual{TagAD1{T},T,1}

# ------------------------------------------------------------------------------------------
# Points 

const FrameFunWrapper = FunctionWrappersWrapper

function _frame_point_fun_wrapper(::Val{O}, ::Type{T}, fun::Function) where {O,T}
    argtypes = (Tuple{T}, Tuple{DualAD1{T}})
    rettypes = (Translation{O,T}, Translation{O,DualAD1{T}})
    return FunctionWrappersWrapper(
        fun, argtypes, rettypes; cache = NoCache(), policy = AllowAll()
    )
end


# ------------------------------------------------------------------------------------------
# Axes

function _frame_axes_fun_wrapper(::Val{O}, ::Type{T}, fun::Function) where {O,T}
    argtypes = (Tuple{T}, Tuple{DualAD1{T}})
    rettypes = (Rotation{O,T}, Rotation{O,DualAD1{T}})
    return FunctionWrappersWrapper(
        fun, argtypes, rettypes; cache = NoCache(), policy = AllowAll()
    )
end
