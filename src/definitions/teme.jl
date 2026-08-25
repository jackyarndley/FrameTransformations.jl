"""
    teme_rot3_gcrf_to_teme(tt_seconds, model=iers2010b)

Compute the direction cosine matrix from GCRF to the True Equator, Mean Equinox
(TEME) frame. `tt_seconds` is J2000 seconds in Terrestrial Time. This follows
Orekit's `TEMEProvider`: TEME is obtained from True-of-Date by a frame rotation
through the full equation of the equinoxes, including complementary terms and
the available EOP nutation-in-longitude correction.

See <https://www.orekit.org/site-orekit-latest/xref/org/orekit/frames/TEMEProvider.html>.

!!! warning
    TEME has no official definition and should be limited to TLE/SGP4 workflows.
"""
function teme_rot3_gcrf_to_teme(tt_seconds::Number, model::IERSModel=iers2010b)
    tt_centuries = tt_seconds / CENTURY2SEC
    correction = eop_δΔψ(model, tt_centuries)
    equation = equation_equinoxes(model, tt_centuries, correction)
    tod_to_teme = angle_to_dcm(equation, :Z)
    return tod_to_teme * iers_rot3_gcrf_to_tod(tt_seconds, model)
end

function _teme_rotation(tt_seconds, model, order::Val)
    rotation = time -> teme_rot3_gcrf_to_teme(time, model)
    return map(DCM, _value_derivatives(rotation, tt_seconds, order))
end

"""Return the GCRF-to-TEME rotation and its first time derivative."""
teme_rot6_gcrf_to_teme(tt_seconds::Number, model::IERSModel=iers2010b) =
    _teme_rotation(tt_seconds, model, Val(2))

"""Return the GCRF-to-TEME rotation and its first two time derivatives."""
teme_rot9_gcrf_to_teme(tt_seconds::Number, model::IERSModel=iers2010b) =
    _teme_rotation(tt_seconds, model, Val(3))

"""Return the GCRF-to-TEME rotation and its first three time derivatives."""
teme_rot12_gcrf_to_teme(tt_seconds::Number, model::IERSModel=iers2010b) =
    _teme_rotation(tt_seconds, model, Val(4))

"""
    add_axes_teme!(frames, name=:TEME, parentid=AXESID_ICRF, id=AXESID_TEME;
        model=iers2010b)

Add True Equator, Mean Equinox axes to `frames`. A prepared conversion from the
frame-system time scale to Terrestrial Time is captured during registration. Missing
kinematic orders are generated with DifferentiationInterface.
"""
function add_axes_teme!(
    frames::FrameSystem,
    name::Symbol=:TEME,
    parentid::Int=AXESID_ICRF,
    id::Int=AXESID_TEME;
    model::IERSModel=iers2010b,
)
    _require_celestial_parent(parentid, "True Equator, Mean Equinox (TEME)")
    to_tt = Tempo.prepare_time_conversion(timescale(frames)(), TT)
    rotation = time -> teme_rot3_gcrf_to_teme(to_tt(time), model)
    return add_axes_rotating!(frames, name, id, parentid, rotation)
end
