
const POINT_CLASSID_ROOT = 0

const POINT_CLASSID_FIXED = 1

const POINT_CLASSID_DYNAMIC = 2

""" 
    add_point!(frames, name, id, axesid, class, funs, parentid=nothing)

Create and add a new point node `name` to `frames` based on the input parameters. 

### Inputs 
- `frames` -- Target frame system 
- `name` -- Point name, must be unique within `frames` 
- `id` -- Point ID, must be unique within `frames`
- `axesid` -- ID of the axes in which the state vector of the point is expressed. 
- `class` -- Point class. 
- `funs` -- `FramePointFunctions` object storing the functions to update the state 
            vectors of the point. It must match the type and order of `frames`
- `parentid` -- NAIF ID of the parent point. Not required only for the root point.

!!! warning 
    This is a low-level function and is NOT meant to be directly used. Instead, to add a point 
    to the frame system, see [`add_point_dynamical!`](@ref), [`add_point_fixedoffset!`](@ref)
    and [`add_point_root!`](@ref).
"""
function add_point!(
    frames::FrameSystem{O, N}, name::Symbol, id::Int, axesid::Int, class::Int,
    funs::FramePointFunctions{O, N}, parentid=nothing
) where {O, N <: Number}

    if has_point(frames, id)
        # Check point with the same id already registered 
        throw(
            ArgumentError(
                "A point with ID $id is already registered in the input frame system.",
            ),
        )
    end

    # Check point with the same name does not already exist 
    if name in map(x -> x.name, points_graph(frames).nodes)
        throw(
            ArgumentError(
                "A point with name=$name is already registed in the input frame system"
            ),
        )
    end

    # Check if the given axes are known in the FrameSystem
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

    # Creates point node 
    pnt = FramePointNode{O, N}(name, class, id, parentid, axesid, funs)

    # Insert new point in the graph
    add_point!(frames, pnt)

    # Connect the new point to the parent point in the graph 
    !isnothing(parentid) && add_edge!(points_graph(frames), parentid, id)

    return nothing
end

"""
    add_point_root!(frames, name, id, axes)

Add root `name` root point with the specified `id` to `frames`.
"""
function add_point_root!(frames::FrameSystem{O, N}, name::Symbol, id::Int, ax) where {O, N}

    # Check for root-point existence 
    if !isempty(points_graph(frames))
        throw(
            ArgumentError("A root-point is already registed in the input frame system.")
        )
    end

    return add_point!(
        frames, name, id, axes_id(frames, ax), POINT_CLASSID_ROOT, FramePointFunctions{O, N}()
    )
end

"""
    add_point_fixedoffset!(frames, name, id, parent, axes, offset::AbstractVector)

Add `point` as a fixed-offset point to `frames`. Fixed points are those whose positions have a 
constant `offset` with respect their `parent` points in the given set of `axes`. Thus, points 
eligible for this class must have null velocity and acceleration with respect to `parent`.
"""
function add_point_fixedoffset!(
    frames::FrameSystem{O, N}, name::Symbol, id::Int, parent, ax,
    offset::AbstractVector{T}
) where {O, N, T}

    if length(offset) != 3
        throw(
            DimensionMismatch(
                "The offset vector should have length 3, but has $(length(offset))."
            ),
        )
    end

    voffset = SVectorNT{3O, N}(SVector(offset...))
    funs = FramePointFunctions{O, N}(t -> voffset, t -> voffset, t -> voffset, t -> voffset)

    return add_point!(
        frames, name, id, axes_id(frames, ax), POINT_CLASSID_FIXED, funs, 
        point_id(frames, parent)
    )
end

""" 
    add_point_dynamical!(frames, name, id, parent, axes, fun, δfun=nothing, δ²fun=nothing, δ³fun=nothing)

Add `point` as a time point to `frames`. The state vector for these points depends only on 
time and is computed through the custom functions provided by the user. 

The input functions must accept only time as argument and their outputs must be as follows: 

- **fun**: return a 3-elements vector: position
- **δfun**: return a 6-elements vector: position and velocity
- **δ²fun**: return a 9-elements vector: position, velocity and acceleration
- **δ³fun**: return a 12-elements vector: position, velocity, acceleration and jerk

If `δfun`, `δ²fun` or `δ³fun` are not provided, they are computed with automatic differentiation. 

!!! warning 
    It is expected that the input functions and their ouputs have the correct signature. This 
    function does not perform any checks on whether the returned vectors have the appropriate 
    dimensions. 
"""
function add_point_dynamical!(
    frames::FrameSystem{O, N}, name::Symbol, id::Int, parent, ax,
    fun, δfun = nothing, δ²fun = nothing, δ³fun = nothing,
) where {O, N}

    for (order, fcn) in enumerate([δfun, δ²fun, δ³fun])
        if (O < order + 1 && !isnothing(fcn))
            @warn "ignoring $fcn, frame system order is less than $(order+1)"
        end
    end

    if O > 1 && isnothing(δfun)
        δfun = t->vcat(fun(t), D¹(fun, t))[1:min(6, 3O)]
    end

    if O > 2 && isnothing(δ²fun)
        δ²fun = t->vcat(δfun(t), D²(fun, t))[1:min(9, 3O)]
    end

    if O > 3 && isnothing(δ³fun)
        δ³fun = t->vcat(δ²fun(t), D³(fun, t))[1:3O]
    end

    tmp = [fun, δfun, δ²fun, δ³fun]

    funs = FramePointFunctions{O, N}(
        [t->vcat(tmp[i](t), zeros(3(O - i))) for i in 1:O]...
    )
    
    return add_point!(
        frames, name, id, axes_id(frames, ax), POINT_CLASSID_DYNAMIC, funs, 
        point_id(frames, parent)
    )
end

"""
    add_point_alias!(frames, target, alias::Symbol)
    add_point_alias!(frames, target, alias::Symbol)

Add a name `alias` to a `target` point registered in `frames`.
"""
function add_point_alias!(frames::FrameSystem{O, N}, target, alias::Symbol) where {O, N}
    if !has_point(frames, target)
        throw(
            ErrorException(
                "no point with ID $target registered in the given frame system"
            )
        )
    end

    if alias in keys(points(frames))
        throw(
            ErrorException(
                "point with name $alias already present in the given frame system"
            )
        )
    end

    push!(points(frames), Pair(alias, point_id(frames, target)))
    nothing
end