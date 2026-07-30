
"""
    add_axes_bcrtod!(fr, name, id, center; deriv=false)

Add Body-Centered Rotating (BCR), True-of-Date (TOD) axes with `name` and `id` to `fr`. 
The center point (i.e., the reference body) is `center`.

These axes are the equivalent of SPICE's `IAU_<BODY_NAME>` frames.

!!! warning 
    The parent axes are automatically set to the ICRF (ID = $(AXESID_ICRF)). If the 
    ICRF is not defined in `fr`, an error is thrown.

### See also 
See also [`add_axes_rotating!`](@ref), [`add_axes_bci2000!`](@ref).
"""
function add_axes_bcrtod!(
    fr::FrameSystem, name::Symbol, id::Int, center; deriv::Bool=false
)
    # create val dispatch
    cid = point_id(fr, center)
    vid = Val(cid)

    if length(methods(body_rotational_elements, [Number, Val{cid}])) != 1
        throw(
            MethodError(
                body_rotational_elements,
                "there must be a unique method defined for $cid"
            )
        )
    end

    if !(has_axes(fr, AXESID_ICRF))
        throw(
            ErrorException(
                "Body-Centered Inertial (BCI) at J2000 axes can only be defined" *
                " w.r.t. the ICRF (ID = $(AXESID_ICRF)), which is not defined in" *
                " the current frames graph."
            )
        )
    end

    if !deriv
        add_axes_rotating!(fr, name, id, AXESID_ICRF, t -> _bcrtod(t, vid))
    else
        add_axes_rotating!(
            fr, name, id, AXESID_ICRF,
            t -> _bcrtod(t, vid),
            t -> _bcrtod_derivative1(t, vid),
            t -> _bcrtod_derivative2(t, vid),
            t -> _bcrtod_derivative3(t, vid),
        )
    end
end

"""
    add_axes_bci2000!(fr, axes::AbstractFrameAxes, center, data)

Add Body-Centered Inertial (BCI) axes at J2000 with `name` and `id` to `fr`. 
The center point (i.e., the reference body) is `center`.

!!! warning 
    The parent axes are automatically set to the ICRF (ID = $(AXESID_ICRF)). If the 
    ICRF is not defined in `fr`, an error is thrown.
    
### See also 
See also [`add_axes_fixedoffset!`](@ref), [`add_axes_bcrtod!`](@ref).
"""
function add_axes_bci2000!(fr::FrameSystem, name::Symbol, id::Int, center)

    # create val dispatch
    cid = point_id(fr, center)
    vid = Val(cid)

    if length(methods(body_rotational_elements, [Number, Val{cid}])) != 1
        throw(
            ErrorException(
                "there must be a unique method defined for $cid"
            )
        )
    end

    if !(has_axes(fr, AXESID_ICRF))
        throw(
            ErrorException(
                "Body-Centered Inertial (BCI) at J2000 axes can only be defined" *
                " w.r.t. the ICRF (ID = $(AXESID_ICRF)), which is not defined in" *
                " the current frames graph."
            )
        )
    end

    # evaluate rotational elements and build the DCM 
    right_ascension_j2000, declination_j2000, _ =
        body_rotational_elements(0, vid)
    dcm = angle_to_dcm(
        pi / 2 + right_ascension_j2000,
        pi / 2 - declination_j2000,
        :ZX,
    )

    # insert the new axes
    return add_axes_fixedoffset!(fr, name, id, AXESID_ICRF, dcm)

end

function _bcrtod(seconds, val)
    right_ascension, declination, prime_meridian =
        body_rotational_elements(seconds / CENTURY2SEC, val)
    return angle_to_dcm(
        pi / 2 + right_ascension,
        pi / 2 - declination,
        prime_meridian,
        :ZXZ,
    )
end

function _bcrtod_derivative1(seconds, val)
    rotation_function = epoch -> _bcrtod(epoch, val)
    return (
        rotation_function(seconds),
        derivative1(rotation_function, seconds),
    )
end

function _bcrtod_derivative2(seconds, val)
    rotation_function = epoch -> _bcrtod(epoch, val)
    return (
        rotation_function(seconds),
        derivative1(rotation_function, seconds),
        derivative2(rotation_function, seconds),
    )
end

function _bcrtod_derivative3(seconds, val)
    rotation_function = epoch -> _bcrtod(epoch, val)
    return (
        rotation_function(seconds),
        derivative1(rotation_function, seconds),
        derivative2(rotation_function, seconds),
        derivative3(rotation_function, seconds),
    )
end
