
const POINT_CLASSID_ROOT = 0

const POINT_CLASSID_GENERIC = -1

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
    if name in map(x -> x.name, get_points(frames).nodes)
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
        if !isempty(get_points(frames))
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
    !isnothing(parentid) && add_edge!(get_points(frames), parentid, id)

    return nothing
end


function add_point_root!(
    frames::FrameSystem{O, N}, name::Symbol, id::Int, axesid::Int
) where {O, N}

    # Check for root-point existence 
    if !isempty(get_points(frames))
        throw(
            ArgumentError("A root-point is already registed in the input frame system.")
        )
    end

    return add_point!(
        frames, name, id, axesid, POINT_CLASSID_ROOT, FramePointFunctions{O, N}(), id
    )
end
