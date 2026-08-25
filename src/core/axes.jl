"""
    add_axes!(frames, name::Symbol, id::Int, class::Int, funs, parentid)

Add a new axes node to `frames`.

### Inputs 
- `frames` -- Target frame system 
- `name` -- Axes name, must be unique within `frames` 
- `id` -- Axes ID, must be unique within `frames`
- `funs` -- `FrameAxesFunctions` object storing the functions to compute the DCM and, 
            eventually, its time derivatives. 
- `parentid` -- Axes ID of the parent axes. Not required only for the root axes.

!!! warning 
    This is a low-level function and is NOT meant to be directly used. Instead, to add a set of
    axes to the frame system, see [`add_axes_projected!`](@ref), [`add_axes_rotating!`](@ref) 
    and [`add_axes_fixedoffset!`](@ref).
"""
function add_axes!(
    frames::FrameSystem{O,T}, name::Symbol, id::Int,
    funs::FrameAxesFunctions{O,T}=FrameAxesFunctions{O,T}(),
    parentid=nothing
) where {O,T<:Number}

    if has_axes(frames, id)
        # Check if a set of axes with the same ID is already registered within 
        # the given frame system 
        throw(
            ArgumentError(
                "Axes with ID $id are already registered in the frame system."
            ),
        )
    end

    if haskey(axes_alias(frames), name)
        # Check if axes with the same name also do not already exist
        throw(
            ArgumentError(
                "Axes with name=$name are already registered in the frame system."
            ),
        )
    end

    if !isnothing(parentid)
        # Check if the root axes is not present
        isempty(axes_graph(frames)) && throw(ArgumentError("Missing root axes."))

        # Check if the parent axes are registered in frame 
        if !has_axes(frames, parentid)
            throw(
                ArgumentError(
                    "The specified parent axes with ID $parentid are not " *
                    "registered in the frame system.",
                ),
            )
        end
    else
        # Check if axes are already present 
        !isempty(axes_graph(frames)) && throw(ArgumentError("Root axes already registed."))

        # Root axes
        parentid = id
    end


    # Create node and insert into the graph
    node = FrameAxesNode{O,T}(name, id, parentid, funs)

    # Insert new point in the graph
    add_axes!(frames, node)

    # Connect the new axes to the parent axes in the graph (skip for root)
    parentid != id && add_edge!(axes_graph(frames), parentid, id)

    return nothing
end

"""
    add_axes_fixedoffset!(frames, name::Symbol, id::Int, parent, dcm:DCM)
   
Add axes `name` with id `id` to `frames` with a fixed-offset from `parent`. 
Fixed offset axes have a constant orientation with respect to their `parent` axes, 
represented by `dcm`, a Direction Cosine Matrix (DCM).

### See also 
See also [`add_axes!`](@ref).
"""
function add_axes_fixedoffset!(
    frames::FrameSystem{O,T}, name::Symbol, id::Int, parent, dcm::DCM{T}
) where {O,T}

    rotation = t -> one(t) * dcm
    funs = FrameAxesFunctions{O,T}(rotation)
    add_axes!(frames, name, id, funs, axes_id(frames, parent))
end

"""
    add_axes_projected!(frames, name, id, parent, fun)

Add inertial axes `name` and id `id` as a set of projected axes to `frames`. The axes relation 
to the `parent` axes are given by a `fun`. 

Projected axes are similar to rotating axis, except that all the positions, velocity, etc ... 
are rotated by the 0-order rotation (i.e. the derivatives of the rotation matrix are null, 
despite the rotation depends on time).

### See also 
See also [`add_axes!`](@ref).
"""
function add_axes_projected!(
    frames::FrameSystem{O,T}, name::Symbol, id::Int, parent, fun
) where {O,T}
    funs = FrameAxesFunctions{O,T}(fun)
    add_axes!(frames, name, id, funs, axes_id(frames, parent))
end

function _next_rotation_function(previous, ::Val{N}) where {N}
    return function (t)
        lower = Rotation{N - 1}(previous(t))
        next = derivative1(epoch -> Rotation{N - 1}(previous(epoch))[N - 1], t)
        return Rotation((lower.m..., DCM(next)))
    end
end

function _rotation_from_base(rotation3, order::Val)
    return t -> Rotation(map(DCM, _value_derivatives(rotation3, t, order)))
end

function _rotation_from_second(rotation6, ::Val{3})
    return _next_rotation_function(t -> Rotation{2}(rotation6(t)), Val(3))
end

function _rotation_from_second(rotation6, ::Val{4})
    return function (t)
        lower = Rotation{2}(rotation6(t))
        last_component = epoch -> Rotation{2}(rotation6(epoch))[2]
        acceleration = derivative1(last_component, t)
        jerk = derivative2(last_component, t)
        return Rotation((lower.m..., DCM(acceleration), DCM(jerk)))
    end
end

function _rotation_functions(::Val{1}, rotation3, rotation6, rotation9, rotation12)
    return (t -> Rotation{1}(rotation3(t)),)
end

function _rotation_functions(::Val{2}, rotation3, rotation6, rotation9, rotation12)
    lower = _rotation_functions(Val(1), rotation3, rotation6, rotation9, rotation12)
    next = isnothing(rotation6) ?
        _rotation_from_base(rotation3, Val(2)) :
        (t -> Rotation{2}(rotation6(t)))
    return (lower..., next)
end

function _rotation_functions(::Val{3}, rotation3, rotation6, rotation9, rotation12)
    lower = _rotation_functions(Val(2), rotation3, rotation6, rotation9, rotation12)
    next = if !isnothing(rotation9)
        t -> Rotation{3}(rotation9(t))
    elseif !isnothing(rotation6)
        _rotation_from_second(rotation6, Val(3))
    else
        _rotation_from_base(rotation3, Val(3))
    end
    return (lower..., next)
end

function _rotation_functions(::Val{4}, rotation3, rotation6, rotation9, rotation12)
    lower = _rotation_functions(Val(3), rotation3, rotation6, rotation9, rotation12)
    next = if !isnothing(rotation12)
        t -> Rotation{4}(rotation12(t))
    elseif !isnothing(rotation9)
        _next_rotation_function(t -> Rotation{3}(rotation9(t)), Val(4))
    elseif !isnothing(rotation6)
        _rotation_from_second(rotation6, Val(4))
    else
        _rotation_from_base(rotation3, Val(4))
    end
    return (lower..., next)
end

"""
    add_axes_rotating!(frames, name::Symbol, id::Int, parent, rotation3;
        rotation6=nothing, rotation9=nothing, rotation12=nothing)

Add a set of rotating axes to `frames`. Each function accepts time and returns the
cumulative rotation through its named order: `rotation3` returns a DCM, while
`rotation6`, `rotation9`, and `rotation12` include one, two, or three time
derivatives. Explicit cumulative-order functions take priority; missing orders are
generated by differentiating the highest available lower-order representation.
"""
function add_axes_rotating!(
    frames::FrameSystem{O,T}, name::Symbol, id::Int, parent, rotation3;
    rotation6=nothing, rotation9=nothing, rotation12=nothing,
) where {O,T}

    for (derivative_order, function_object) in enumerate(
        (rotation6, rotation9, rotation12)
    )
        if O < derivative_order + 1 && !isnothing(function_object)
            @warn "ignoring $function_object, frame system order is less than $(derivative_order + 1)"
        end
    end

    functions = _rotation_functions(
        Val(O), rotation3, rotation6, rotation9, rotation12)
    funs = FrameAxesFunctions{O,T}(functions...)

    return add_axes!(frames, name, id, funs, axes_id(frames, parent))
end

"""
    add_axes_alias!(frames, target, alias::Symbol)

Add a name `alias` to a `target` axes registered in `frames`.
"""
function add_axes_alias!(frames::FrameSystem{O,T}, target, alias::Symbol) where {O,T}
    if !has_axes(frames, target)
        throw(
            ErrorException(
                "no axes with ID $target registered in the given frame system"
            )
        )
    end

    if haskey(axes_alias(frames), alias)
        throw(
            ErrorException(
                "axes with name $alias already present in the given frame system"
            )
        )
    end

    push!(axes_alias(frames), Pair(alias, axes_id(frames, target)))
    nothing
end
