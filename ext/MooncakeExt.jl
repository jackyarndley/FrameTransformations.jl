module MooncakeExt

using FrameTransformations: FrameSystem,
    PreparedRotation, PreparedTranslation, PreparedDirection,
    CompiledRotation, CompiledTranslation, CompiledDirection,
    vector3, vector6, vector9, vector12,
    rotation3, rotation6, rotation9, rotation12,
    direction3, direction6, direction9, direction12

using FunctionWrappers: FunctionWrapper
using FunctionWrappersWrappers: FunctionWrappersWrapper

import Mooncake
using Mooncake: @from_chainrules, DefaultCtx

# ==========================================================================================
# Declare FrameSystem and its internal opaque types as non-differentiable
# ==========================================================================================

# FunctionWrapper and FunctionWrappersWrapper are compiled C function pointers — opaque
# to reverse-mode AD. Declaring them as NoTangent prevents Mooncake from trying to
# construct tangent types for them (which would fail).
Mooncake.tangent_type(::Type{<:FunctionWrapper}) = Mooncake.NoTangent
Mooncake.tangent_type(::Type{<:FunctionWrappersWrapper}) = Mooncake.NoTangent

# FrameSystem is a structural container (graph of axes/points with function wrappers).
# It is never differentiated through — only `t` is differentiable.
Mooncake.tangent_type(::Type{<:FrameSystem}) = Mooncake.NoTangent
const PreparedOrCompiled = Union{
    PreparedRotation,
    PreparedTranslation,
    PreparedDirection,
    CompiledRotation,
    CompiledTranslation,
    CompiledDirection,
}
Mooncake.tangent_type(::Type{<:PreparedOrCompiled}) = Mooncake.NoTangent

# ==========================================================================================
# Register ChainRulesCore rrules with Mooncake
# ==========================================================================================

# Vector functions: vector3, vector6, vector9, vector12
@from_chainrules DefaultCtx Tuple{typeof(vector3), FrameSystem, Any, Any, Any, T} where {T<:Number}
@from_chainrules DefaultCtx Tuple{typeof(vector6), FrameSystem, Any, Any, Any, T} where {T<:Number}
@from_chainrules DefaultCtx Tuple{typeof(vector9), FrameSystem, Any, Any, Any, T} where {T<:Number}
@from_chainrules DefaultCtx Tuple{typeof(vector12), FrameSystem, Any, Any, Any, T} where {T<:Number}

# Rotation functions: rotation3, rotation6, rotation9, rotation12
@from_chainrules DefaultCtx Tuple{typeof(rotation3), FrameSystem, Any, Any, T} where {T<:Number}
@from_chainrules DefaultCtx Tuple{typeof(rotation6), FrameSystem, Any, Any, T} where {T<:Number}
@from_chainrules DefaultCtx Tuple{typeof(rotation9), FrameSystem, Any, Any, T} where {T<:Number}
@from_chainrules DefaultCtx Tuple{typeof(rotation12), FrameSystem, Any, Any, T} where {T<:Number}

# Direction functions: direction3, direction6, direction9, direction12
@from_chainrules DefaultCtx Tuple{typeof(direction3), FrameSystem, Symbol, Any, T} where {T<:Number}
@from_chainrules DefaultCtx Tuple{typeof(direction6), FrameSystem, Symbol, Any, T} where {T<:Number}
@from_chainrules DefaultCtx Tuple{typeof(direction9), FrameSystem, Symbol, Any, T} where {T<:Number}
@from_chainrules DefaultCtx Tuple{typeof(direction12), FrameSystem, Symbol, Any, T} where {T<:Number}

# Prepared and compiled callables differentiate only with respect to scalar time.
@from_chainrules DefaultCtx Tuple{C, T} where {C<:PreparedRotation,T<:Number}
@from_chainrules DefaultCtx Tuple{C, T} where {C<:PreparedTranslation,T<:Number}
@from_chainrules DefaultCtx Tuple{C, T} where {C<:PreparedDirection,T<:Number}
@from_chainrules DefaultCtx Tuple{C, T} where {C<:CompiledRotation,T<:Number}
@from_chainrules DefaultCtx Tuple{C, T} where {C<:CompiledTranslation,T<:Number}
@from_chainrules DefaultCtx Tuple{C, T} where {C<:CompiledDirection,T<:Number}

end # module
