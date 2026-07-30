
"""
    add_axes_topocentric!(frames, name::Symbol, id::Int, parent, longitude, latitude, mount)

Add topocentric axes to `frames` at a specified location and mounting.

The orientation relative to the parent axes `parent` is defined through `longitude`,
geodetic `latitude`, and the mounting type `mount`, which may be any of the following:

- `:NED` (North, East, Down): the X-axis points North, the Y-axis is directed eastward and 
    the Z-axis points inwards towards the nadir.
- `:SEZ` (South, East, Zenith): the X-axis points South, the Y-axis is directed East, and 
    the Z-axis points outwards towards the zenith.
- `:ENU` (East, North, Up): the X-axis points East, the Y-axis is directed North and the 
    Z-axis points outwards towards the zenith. 

!!! warning 
    The parent axes must be a set of body-fixed reference axes. This is under user 
    resposibility. 
"""
function add_axes_topocentric!(
    frames::FrameSystem, name::Symbol, id::Int, parent,
    longitude::Number, latitude::Number, mount::Symbol
)
    if mount == :NED
        dcm = angle_to_dcm(longitude, -latitude - pi / 2, :ZY)
    elseif mount == :SEZ
        dcm = angle_to_dcm(longitude, pi / 2 - latitude, :ZY)
    elseif mount == :ENU
        dcm = angle_to_dcm(longitude + pi / 2, pi / 2 - latitude, :ZX)
    else
        throw(ArgumentError("$mount is not a supported topocentric mounting type."))
    end
    pid = axes_id(frames, parent)
    return add_axes_fixedoffset!(frames, name, id, pid, dcm)
end

@fastmath function _geodetic_to_position(
    height::Number,
    longitude::Number,
    latitude::Number,
    radius::Number,
    flattening::Number,
)
    # Get eccentricity from flattening 
    eccentricity_squared = (2 - flattening) * flattening

    sin_latitude, cos_latitude = sincos(latitude)
    sin_longitude, cos_longitude = sincos(longitude)

    prime_vertical_radius =
        radius / sqrt(1 - eccentricity_squared * sin_latitude^2)
    radial = (prime_vertical_radius + height) * cos_latitude
    polar = (1 - eccentricity_squared) * prime_vertical_radius

    return SA[
        radial * cos_longitude,
        radial * sin_longitude,
        (polar + height) * sin_latitude,
    ]
end

"""
    add_point_surface!(frames, name::Symbol, pointid::Int, parent, axes, 
        longitude::Number, latitude::Number, radius::Number,
        flattening::Number=0.0, height::Number=0.0)

Add `point` to `frames` as a fixed point on the surface of the `parent` point body. 
The relative position is specified by `longitude`, geodetic `latitude`, the reference
`radius` of the ellipsoid, and its `flattening`. `height` defaults to zero.

!!! warning 
    Axes used here must be a set of body-fixed reference axes for the body represented by 
    `parentid`. This is under user resposibility. 
"""
function add_point_surface!(
    frames::FrameSystem, name::Symbol, id::Int, parent, axes,
    longitude::Number,
    latitude::Number,
    radius::Number,
    flattening::Number=0.0,
    height::Number=0.0,
)
    pos = _geodetic_to_position(
        height, longitude, latitude, radius, flattening
    )
    pid = point_id(frames, parent)
    axid = axes_id(frames, axes)
    return add_point_fixedoffset!(frames, name, id, pid, axid, pos)
end
