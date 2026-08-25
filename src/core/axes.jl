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

    empty!(frames._axes_nodes)
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

    rotation = t -> dcm
    funs = _ordered_frame_axes_functions(
        Val(O), T, ntuple(_ -> rotation, Val(O)))
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
    funs = _ordered_frame_axes_functions(
        Val(O), T, ntuple(_ -> fun, Val(O)))
    add_axes!(frames, name, id, funs, axes_id(frames, parent))
end

"""
    add_axes_rotating!(frames, name::Symbol, id::Int, parent, fun,
        first_derivative=nothing, second_derivative=nothing, third_derivative=nothing)
   
Add `axes` as a set of rotating axes to `frames`. The orientation of these axes depends only 
on time and is computed through the custom functions provided by the user. 

The input functions must accept only time as argument and their outputs must be as follows: 

- `fun`: return a Direction Cosine Matrix (DCM).
- `first_derivative`: return the DCM and its first time derivative.
- `second_derivative`: return the DCM and its first two time derivatives.
- `third_derivative`: return the DCM and its first three time derivatives.

Missing derivative functions are computed via automatic differentiation.

!!! warning 
    It is expected that the input functions and their outputs have the correct signature. This 
    function does not perform any checks on the output types. 
"""
function add_axes_rotating!(
    frames::FrameSystem{O,T}, name::Symbol, id::Int, parent, fun,
    first_derivative=nothing, second_derivative=nothing, third_derivative=nothing,
) where {O,T}

    for (derivative_order, function_object) in enumerate(
        (first_derivative, second_derivative, third_derivative)
    )
        if O < derivative_order + 1 && !isnothing(function_object)
            @warn "ignoring $function_object, frame system order is less than $(derivative_order + 1)"
        end
    end

    functions = (
        t -> Rotation{1}(fun(t)),

        # First derivative 
        if isnothing(first_derivative)
            t -> Rotation{2}(fun(t), derivative1(fun, t))
        else
            t -> Rotation{2}(first_derivative(t))
        end,

        # Second derivative 
        if isnothing(second_derivative)
            (
                if isnothing(first_derivative)
                    t -> Rotation{3}(
                        fun(t), derivative1(fun, t), derivative2(fun, t)
                    )
                else
                    t -> Rotation{3}(
                        first_derivative(t)..., derivative2(fun, t)
                    )
                end
            )
        else
            t -> Rotation{3}(second_derivative(t))
        end,

        # Third derivative 
        if isnothing(third_derivative)
            (
                if isnothing(second_derivative)
                    (
                        if isnothing(first_derivative)
                            t -> Rotation{4}(
                                fun(t),
                                derivative1(fun, t),
                                derivative2(fun, t),
                                derivative3(fun, t),
                            )
                        else
                            t -> Rotation{4}(
                                first_derivative(t)...,
                                derivative2(first_derivative, t)...,
                            )
                        end
                    )
                else
                    t -> Rotation{4}(
                        second_derivative(t)..., derivative3(fun, t)
                    )
                end
            )
        else
            t -> Rotation{4}(third_derivative(t))
        end,
    )

    funs = _ordered_frame_axes_functions(Val(O), T, functions)

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
