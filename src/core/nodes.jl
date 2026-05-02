# ------------------------------------------------------------------------------------------
# POINTS
# ------------------------------------------------------------------------------------------

# ------
# Functions

struct FramePointFunctions{O,T,FW<:FrameFunWrapper}
    fun::NTuple{O,FW}
end

Base.getindex(pf::FramePointFunctions, i) = pf.fun[i]
@inline Base.getindex(pf::FramePointFunctions, ::Val{I}) where {I} = getfield(pf.fun, I)

function FramePointFunctions{T}(funs::Vararg{Function,O}) where {T,O}
    wrappers = ntuple(i -> _frame_point_fun_wrapper(Val(O), T, funs[i]), Val(O))
    return FramePointFunctions{O,T,eltype(typeof(wrappers))}(wrappers)
end

function FramePointFunctions{O,T}(funs::Function...) where {O,T}
    O > length(funs) && throw(ArgumentError("required at least $O functions."))
    wrappers = ntuple(i -> _frame_point_fun_wrapper(Val(O), T, funs[i]), Val(O))
    return FramePointFunctions{O,T,eltype(typeof(wrappers))}(wrappers)
end

function FramePointFunctions{O,T}(fun::Function) where {O,T}
    wrappers = ntuple(_ -> _frame_point_fun_wrapper(Val(O), T, fun), Val(O))
    return FramePointFunctions{O,T,eltype(typeof(wrappers))}(wrappers)
end

function FramePointFunctions{O,T}() where {O,T}
    return FramePointFunctions{O,T}(t -> Translation{O,T}())
end


# ------
# Node

"""
    FramePointNode{O, T} <: AbstractJSMDGraphNode

Define a frame system point.

### Fields
- `name` -- point name 
- `id` -- ID of the point
- `parentid` -- ID of the parent point 
- `axesid` -- ID of the axes in which the point coordinates are expressed 
- `f` -- `FramePointFunctions` container 
"""
struct FramePointNode{O,T<:Number} <: AbstractJSMDGraphNode
    name::Symbol
    id::Int
    parentid::Int
    axesid::Int

    # internals
    f::FramePointFunctions{O,T}
end

get_node_id(p::FramePointNode{O,T}) where {O,T} = p.id

function Base.show(io::IO, p::FramePointNode{O,T}) where {O,T}
    pstr = "FramePointNode{$O, $T}(name=$(p.name)"
    pstr *= ", id=$(p.id), axesid=$(p.axesid)"
    p.parentid == p.id || (pstr *= ", parent=$(p.parentid)")
    pstr *= ")"
    return print(io, pstr)
end

const PointsGraph{O,T} = MappedNodeGraph{FramePointNode{O,T},SimpleGraph{Int}}

# ------------------------------------------------------------------------------------------
# AXES
# ------------------------------------------------------------------------------------------

# ------
# Functions

struct FrameAxesFunctions{O,T,FW<:FrameFunWrapper}
    fun::NTuple{O,FW}
end

Base.getindex(pf::FrameAxesFunctions, i) = pf.fun[i]
@inline Base.getindex(pf::FrameAxesFunctions, ::Val{I}) where {I} = getfield(pf.fun, I)

function FrameAxesFunctions{T}(funs::Vararg{Function,O}) where {T,O}
    wrappers = ntuple(i -> _frame_axes_fun_wrapper(Val(O), T, funs[i]), Val(O))
    return FrameAxesFunctions{O,T,eltype(typeof(wrappers))}(wrappers)
end

function FrameAxesFunctions{O,T}(funs::Function...) where {O,T}
    O > length(funs) && throw(ArgumentError("required at least $O functions."))
    wrappers = ntuple(i -> _frame_axes_fun_wrapper(Val(O), T, funs[i]), Val(O))
    return FrameAxesFunctions{O,T,eltype(typeof(wrappers))}(wrappers)
end

function FrameAxesFunctions{O,T}(fun::Function) where {O,T}
    wrappers = ntuple(_ -> _frame_axes_fun_wrapper(Val(O), T, fun), Val(O))
    return FrameAxesFunctions{O,T,eltype(typeof(wrappers))}(wrappers)
end

function FrameAxesFunctions{O,T}() where {O,T}
    return FrameAxesFunctions{O,T}(t -> Rotation{O,T}(one(T)I))
end

# ------
# Node

"""
    FrameAxesNode{O, T} <: AbstractJSMDGraphNode

Define a set of axes.

### Fields
- `name` -- axes name 
- `id` -- axes ID (equivalent of NAIFId for axes)
- `parentid` -- ID of the parent axes 
- `f` -- `FrameAxesFunctions` container 
"""
struct FrameAxesNode{O,T<:Number} <: AbstractJSMDGraphNode
    name::Symbol
    id::Int
    parentid::Int

    # internals
    f::FrameAxesFunctions{O,T}
end

get_node_id(ax::FrameAxesNode{O,T}) where {O,T} = ax.id

function Base.show(io::IO, ax::FrameAxesNode{O,T}) where {O,T}
    pstr = "FrameAxesNode{$O, $T}(name=$(ax.name), id=$(ax.id)"
    ax.parentid == ax.id || (pstr *= ", parent=$(ax.parentid)")
    pstr *= ")"
    return print(io, pstr)
end

const AxesGraph{O,T} = MappedNodeGraph{FrameAxesNode{O,T},SimpleGraph{Int}}

# ------------------------------------------------------------------------------------------
# DIRECTIONS
# ------------------------------------------------------------------------------------------

const DirectionFunctions{O,T} = FramePointFunctions{O,T}

"""
    DirectionDefinition{O, T}

Define a new direction.

### Fields
- `name` -- direction name 
- `id` -- direction ID
- `f` -- `DirectionFunctions` container 
"""
struct DirectionDefinition{O,T}
    name::Symbol
    id::Int
    axesid::Int

    # internals
    f::DirectionFunctions{O,T}
end

function Base.show(io::IO, d::DirectionDefinition{O,T}) where {O,T}
    return print(io, "DirectionDefinition{$O, $T}(name=$(d.name), id=$(d.id), axesid=$(d.axesid))")
end
