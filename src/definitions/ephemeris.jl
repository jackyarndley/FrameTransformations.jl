
"""
    add_point_ephemeris!(fr::FrameSystem{O, T}, eph::AbstractEphemerisProvider, 
        name::Symbol, id::Int) where {O, T} 

Add a point coming from an `AbstractEphemerisProvider` subtype. 
The point is identifies by the `id` and a have a user defined `name`.

!!! warning 
    This is an interface only, concrete subtypes of `AbstractEphemerisProvider` requires 
    an proper implementation.
"""
@interface function add_point_ephemeris!(
    ::FrameSystem{O,T}, ::AbstractEphemerisProvider, ::Symbol, ::Int
) where {O,T} end

"""
    add_axes_ephemeris!(frames::FrameSystem{O, T}, eph::AbstractEphemerisProvider, 
        name::Symbol, id::Int, seq::Symbol) where {O, T} 

Add axes coming from an `AbstractEphemerisProvider` subtype to `frames`. 
The axes are identifies by the `id` and a have a user defined `name`. The rotation matrix 
is build using the rotation sequence specified in `seq`. The axes type is specified by `class`.

!!! warning 
    This is an interface only, concrete subtypes of `AbstractEphemerisProvider` requires 
    an proper implementation.
"""
@interface function add_axes_ephemeris!(
    ::FrameSystem{O,T}, ::AbstractEphemerisProvider, ::Symbol, ::Int, ::Symbol
) where {O,T} end

"""
    add_point_ephemeris!(fr, eph::AbstractEphemerisProvider, book::Dict{Int, Symbol})

Add all points found in the `eph` and using id-names relationships specified in `book`.
"""
function add_point_ephemeris!(
    fr::FrameSystem{O,T}, eph::AbstractEphemerisProvider, book::Dict{Int,Symbol}
) where {O,T}
    records = ephem_position_records(eph)
    for id in sort(ephem_available_points(eph))
        !haskey(book, id) && @warn "Cannot find point with ID $id in the names book" continue

        if id == 0
            if length(points_graph(fr).nodes) == 0
                # No point registered in the frame system 
                add_point!(fr, book[id], id, AXESID_ICRF)
            end
        else
            # Assume that the SSB is registered in the frame system
            it = findfirst(rec -> rec.target == id, records)
            rec = records[it]
            add_point_ephemeris!(fr, eph, book[id], id)
        end
    end
    nothing
end

function check_point_ephemeris(
    fr::FrameSystem{O,T}, eph::AbstractEphemerisProvider, id::Int
) where {O,T}

    # Check that the kernels contain the ephemeris data for the given naifid
    if !(id in ephem_available_points(eph))
                throw(
                    ErrorException(
                        "Ephemeris data for ID $id is not available in the kernels.",
                    ),
                )
    end

    # Retrieve the ephemerides position records (i.e., the segment descriptors)
    pos_records = ephem_position_records(eph)

    # Retrieve the parent point and the point axes from the ephemeris data 
    parentid, parent_found = 0, false
    axesid, axes_found = 0, false

    for pr in pos_records
        if pr.target == id

            # Update the parent point ID
            if !parent_found
                parentid = pr.center
                parent_found = true

            elseif parentid != pr.center
                throw(
                    ErrorException(
                        "UnambiguityError: at least two set of data with different" *
                        " centers are available for point with ID $id.",
                    ),
                )

            end

            # Update the reference axes ID
            if !axes_found
                axesid = pr.axes
                axes_found = true

            elseif axesid != pr.axes
                throw(
                    ErrorException(
                        "UnambiguityError: at least two set of data with different axes" *
                        " are available for point with ID $id.",
                    ),
                )
            end

        end
    end

    # Check that the default parent is available in the FrameSystem
    if !has_point(fr, parentid)
        throw(
            ErrorException(
                "Ephemeris data for point with ID $id is available with respect" *
                " to point with ID $parentid, which has not yet been defined" *
                " in the input frame system.",
            ),
        )
    end

    # The parent point is automatically inferred from the ephemeris kernels so it will 
    # always have available ephemeris data

    # Checks if the axes are known to the frame system. 
    # This check is also performed by add_point!, but it is reported here because 
    # it provides more specific information for ephemeris points 
    if !has_axes(fr, axesid)
        throw(
            ErrorException(
                "Ephemeris data for point with ID $id is expressed in a set" *
                " of axes with ID $axesid, which are yet to be defined in the" *
                " input frame system.",
            ),
        )
    end

    return parentid, axesid
end

function check_axes_ephemeris(
    frames::FrameSystem{O,N}, eph::AbstractEphemerisProvider, id::Int
) where {O,N}

    # Check that the kernels contain orientation data for the given axes ID 
    if !(id in ephem_available_axes(eph))
        throw(
            ArgumentError(
                "Orientation data for AXESID $id is not available " *
                "in the kernels loaded in the given FrameSystem.",
            ),
        )
    end

    # Retrieve the parent from the ephemeris data 
    orient_records = ephem_orient_records(eph)
    parentid = nothing

    for or in orient_records
        if or.target == id
            if isnothing(parentid)
                parentid = or.axes
            elseif parentid != or.axes
                throw(
                    ErrorException(
                        "UnambiguityError: at least two set of orientation data " *
                        "with different centers are available for axes with ID $id.",
                    ),
                )
            end
        end
    end

    # Check that the default parent is available in the FrameSystem! 
    if !has_axes(frames, parentid)
        throw(
            ArgumentError(
                "Orientation data for axes with ID $id is available " *
                "with respect to axes with ID $parentid, which has not yet been defined " *
                "in the given FrameSystem.",
            ),
        )
    end

    return parentid
end

# Functions for an easier definition and handling of ephemeris axes.
@inline function _axis_generator(::Type{T}, ::Val{:X}) where {T}
    return DCM(SMatrix{3,3,T,9}(
        zero(T), zero(T), zero(T),
        zero(T), zero(T), -one(T),
        zero(T), one(T), zero(T),
    ))
end

@inline function _axis_generator(::Type{T}, ::Val{:Y}) where {T}
    return DCM(SMatrix{3,3,T,9}(
        zero(T), zero(T), one(T),
        zero(T), zero(T), zero(T),
        -one(T), zero(T), zero(T),
    ))
end

@inline function _axis_generator(::Type{T}, ::Val{:Z}) where {T}
    return DCM(SMatrix{3,3,T,9}(
        zero(T), -one(T), zero(T),
        one(T), zero(T), zero(T),
        zero(T), zero(T), zero(T),
    ))
end

@inline function _elementary_rotation(
    angle,
    rate,
    acceleration,
    jerk,
    axis::Val{Axis},
    ::Val{O},
) where {Axis,O}
    matrix = angle_to_dcm(angle, Axis)
    generator = _axis_generator(eltype(matrix), axis)
    generator2 = generator * generator
    first = rate * generator * matrix

    O == 1 && return Rotation(matrix)

    second =
        acceleration * generator * matrix +
        rate^2 * generator2 * matrix
    O == 2 && return Rotation(matrix, first)

    generator3 = generator2 * generator
    third =
        jerk * generator * matrix +
        3 * rate * acceleration * generator2 * matrix +
        rate^3 * generator3 * matrix
    O == 3 && return Rotation(matrix, first, second)

    return Rotation(matrix, first, second, third)
end

@generated function _angles_to_rotation(
    angles,
    ::Val{Sequence},
    ::Val{O},
) where {Sequence,O}
    axes = Symbol.(collect(String(Sequence)))
    length(axes) == 3 || throw(
        ArgumentError("expected a three-axis rotation sequence, got $Sequence")
    )

    rotations = map(1:3) do index
        rate = O >= 2 ? :(angles[$(3 + index)]) : :(zero(angles[$index]))
        acceleration =
            O >= 3 ? :(angles[$(6 + index)]) : :(zero(angles[$index]))
        jerk = O >= 4 ? :(angles[$(9 + index)]) : :(zero(angles[$index]))
        return :(
            _elementary_rotation(
                angles[$index],
                $rate,
                $acceleration,
                $jerk,
                Val($(QuoteNode(axes[index]))),
                Val($O),
            )
        )
    end

    return quote
        @inbounds $(rotations[3]) * $(rotations[2]) * $(rotations[1])
    end
end

@inline angles_to_rot3(angles, sequence::Val) =
    _angles_to_rotation(angles, sequence, Val(1))
@inline angles_to_rot6(angles, sequence::Val) =
    _angles_to_rotation(angles, sequence, Val(2))
@inline angles_to_rot9(angles, sequence::Val) =
    _angles_to_rotation(angles, sequence, Val(3))
@inline angles_to_rot12(angles, sequence::Val) =
    _angles_to_rotation(angles, sequence, Val(4))
