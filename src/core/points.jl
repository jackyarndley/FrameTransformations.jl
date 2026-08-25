""" 
    add_point!(frames, name, id, axesid, class, funs, parentid=nothing)

Create and add a new point node `name` to `frames` based on the input parameters. 

### Inputs 
- `frames` -- Target frame system 
- `name` -- Point name, must be unique within `frames` 
- `id` -- Point ID, must be unique within `frames`
- `axes` -- ID/Name of the axes in which the state vector of the point is expressed. 
- `funs` -- `FramePointFunctions` object storing the functions to update the state 
            vectors of the point.
- `parentid` -- NAIF ID of the parent point. Not required only for the root point.

!!! warning 
    This is a low-level function and is NOT meant to be directly used. Instead, to add a point 
    to the frame system, see [`add_point_dynamical!`](@ref) and [`add_point_fixedoffset!`](@ref).
"""
function add_point!(
    frames::FrameSystem{O,T}, name::Symbol, id::Int, axes,
    funs::FramePointFunctions{O,T}=FramePointFunctions{O,T}(), parentid=nothing
) where {O,T<:Number}

    if has_point(frames, id)
        # Check point with the same id already registered 
        throw(
            ArgumentError(
                "A point with ID $id is already registered in the input frame system.",
            ),
        )
    end

    # Check point with the same name does not already exist 
    if haskey(points_alias(frames), name)
        throw(
            ArgumentError(
                "A point with name=$name is already registed in the input frame system"
            ),
        )
    end

    # Check if the given axes are known in the FrameSystem
    axesid = axes_id(frames, axes)
    if !has_axes(frames, axesid)
        throw(
            ArgumentError(
                "Axes with ID $axesid are not registered in the input frame system"
            ),
        )
    end

    if isnothing(parentid)
        # If a root-point exists, check that a parent has been specified 
        if !isempty(points_graph(frames))
            throw(
                ArgumentError(
                    "A parent point is required because the input frame system " *
                    "already contains a root-point.",
                ),
            )
        end

        parentid = id # Root-point has parentid = id

    else
        # Check that the parent point is registered in frames 
        if !has_point(frames, parentid)
            throw(
                ArgumentError(
                    "The specified parent point with id $parentid is not " *
                    "registered in the input frame system.",
                ),
            )
        end
    end

    # Create node and insert into the graph
    pnt = FramePointNode{O,T}(name, id, parentid, axesid, funs)

    # Insert new point in the graph
    add_point!(frames, pnt)

    # Connect the new point to the parent point in the graph (skip for root)
    parentid != id && add_edge!(points_graph(frames), parentid, id)

    return nothing
end

"""
    add_point_fixedoffset!(frames, name, id, parent, axes, offset::AbstractVector)

Add `point` as a fixed-offset point to `frames`. 

Fixed points are those whose positions have a constant `offset` with respect their `parent` 
points in the given set of `axes`. Thus, points eligible for this class must have null 
velocity and acceleration with respect to `parent`.
"""
function add_point_fixedoffset!(
    frames::FrameSystem{O,T}, name::Symbol, id::Int, parent, ax,
    offset::AbstractVector{N}
) where {O,N,T}

    if length(offset) != 3
        throw(
            DimensionMismatch(
                "The offset vector should have length 3, but has $(length(offset))."
            ),
        )
    end

    fixed_position = SVector{3,T}(offset)
    funs = FramePointFunctions{O,T}(t -> one(t) * fixed_position)

    return add_point!(
        frames, name, id, axes_id(frames, ax), funs, point_id(frames, parent)
    )
end

function _next_translation_function(previous, ::Val{N}) where {N}
    return function (t)
        lower = Translation{N - 1}(previous(t))
        next = derivative1(
            epoch -> Translation{N - 1}(previous(epoch))[N - 1], t)
        return Translation((lower.v..., SVector{3}(next)))
    end
end

function _translation_from_base(position, ::Val{2})
    return t -> _translation_from_parts(
        Val(2), position(t), derivative1(position, t))
end

function _translation_from_base(position, ::Val{3})
    return t -> _translation_from_parts(
        Val(3), position(t), derivative1(position, t), derivative2(position, t))
end

function _translation_from_base(position, ::Val{4})
    return t -> _translation_from_parts(
        Val(4),
        position(t),
        derivative1(position, t),
        derivative2(position, t),
        derivative3(position, t),
    )
end

function _translation_from_second(state6, ::Val{3})
    return _next_translation_function(
        t -> Translation{2}(state6(t)), Val(3))
end

function _translation_from_second(state6, ::Val{4})
    return function (t)
        lower = Translation{2}(state6(t))
        last_component = epoch -> Translation{2}(state6(epoch))[2]
        acceleration = SVector{3}(derivative1(last_component, t))
        jerk = SVector{3}(derivative2(last_component, t))
        return Translation((lower.v..., acceleration, jerk))
    end
end

function _translation_functions(::Val{1}, position, state6, state9, state12)
    return (t -> Translation{1}(position(t)),)
end

function _translation_functions(::Val{2}, position, state6, state9, state12)
    lower = _translation_functions(Val(1), position, state6, state9, state12)
    next = isnothing(state6) ?
        _translation_from_base(position, Val(2)) :
        (t -> Translation{2}(state6(t)))
    return (lower..., next)
end

function _translation_functions(::Val{3}, position, state6, state9, state12)
    lower = _translation_functions(Val(2), position, state6, state9, state12)
    next = if !isnothing(state9)
        t -> Translation{3}(state9(t))
    elseif !isnothing(state6)
        _translation_from_second(state6, Val(3))
    else
        _translation_from_base(position, Val(3))
    end
    return (lower..., next)
end

function _translation_functions(::Val{4}, position, state6, state9, state12)
    lower = _translation_functions(Val(3), position, state6, state9, state12)
    next = if !isnothing(state12)
        t -> Translation{4}(state12(t))
    elseif !isnothing(state9)
        _next_translation_function(t -> Translation{3}(state9(t)), Val(4))
    elseif !isnothing(state6)
        _translation_from_second(state6, Val(4))
    else
        _translation_from_base(position, Val(4))
    end
    return (lower..., next)
end

"""
    add_point_dynamical!(frames, name, id, parent, axes, position;
        state6=nothing, state9=nothing, state12=nothing)

Add a dynamical point to `frames`. Each function accepts time and returns the
cumulative translation through its named order: `position` returns three elements,
while `state6`, `state9`, and `state12` include velocity, acceleration, or jerk.
Explicit cumulative-order functions take priority; missing orders are generated by
differentiating the highest available lower-order representation.
"""
function add_point_dynamical!(
    frames::FrameSystem{O,T}, name::Symbol, id::Int, parent, ax, position;
    state6=nothing, state9=nothing, state12=nothing,
) where {O,T}

    for (derivative_order, function_object) in enumerate(
        (state6, state9, state12)
    )
        if O < derivative_order + 1 && !isnothing(function_object)
            @warn "ignoring $function_object, frame system order is less than $(derivative_order + 1)"
        end
    end

    functions = _translation_functions(
        Val(O), position, state6, state9, state12)
    funs = FramePointFunctions{O,T}(functions...)

    return add_point!(
        frames, name, id, axes_id(frames, ax), funs, point_id(frames, parent)
    )
end

"""
    add_point_alias!(frames, target, alias::Symbol)

Add a name `alias` to a `target` point registered in `frames`.
"""
function add_point_alias!(frames::FrameSystem{O,N}, target, alias::Symbol) where {O,N}
    if !has_point(frames, target)
        throw(
            ErrorException(
                "no point with ID $target registered in the given frame system"
            )
        )
    end

    if haskey(points_alias(frames), alias)
        throw(
            ErrorException(
                "point with name $alias already present in the given frame system"
            )
        )
    end

    push!(points_alias(frames), Pair(alias, point_id(frames, target)))
    nothing
end
