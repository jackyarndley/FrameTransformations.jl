module FrameTransformations

using LinearAlgebra
using StaticArrays
using ReferenceFrameRotations
using FunctionWrappersWrappers: FunctionWrappersWrapper, AllowAll, NoCache
import DifferentiationInterface as DI
import ForwardDiff # Activates DifferentiationInterface's ForwardDiff extension.

using JSMDUtils.Math:
    arcsec2rad,
    unitvec,
    cross3,
    cross6,
    cross9,
    cross12

using JSMDInterfaces.Graph: AbstractJSMDGraphNode,
       add_edge!, add_vertex!, get_path, has_vertex

using SMDGraphs: MappedNodeGraph, SimpleGraph, MappedGraph,
       get_mappedid, get_mappednode, get_node, get_path

import SMDGraphs: get_node_id

import Tempo
using Tempo: AbstractTimeScale, Epoch, j2000s, BarycentricDynamicalTime,
       TT, CENTURY2SEC

using JSMDInterfaces.Ephemeris: AbstractEphemerisProvider,
       ephem_position_records, ephem_available_points,
       ephem_orient_records, ephem_available_axes

using JSMDInterfaces.Interface: @interface

using JSMDInterfaces.Bodies: body_rotational_elements

using IERSConventions: iers_bias, iers_obliquity, equation_equinoxes, eop_δΔψ,
       iers_rot3_gcrf_to_itrf, iers_rot6_gcrf_to_itrf,
       iers_rot9_gcrf_to_itrf, iers_rot12_gcrf_to_itrf,
       iers_rot3_gcrf_to_mod, iers_rot3_gcrf_to_tod, iers_rot3_gcrf_to_gtod,
       iers_rot3_gcrf_to_pef, iers_rot3_gcrf_to_cirf, iers_rot3_gcrf_to_tirf,
       IERSModel, iers2010a, iers2010b, iers1996

# ==========================================================================================
# Core
# ==========================================================================================

# Low-level types and aliases
export Translation, Rotation

include("core/translation.jl")
include("core/rotation.jl")
include("core/ad.jl")

# Frame system 
export FrameSystem,
       order, timescale, points_graph, axes_graph, points_alias, axes_alias, directions,
       has_axes, has_point, has_direction,
       point_id, axes_id

include("core/nodes.jl")
include("core/graph.jl")

# Helper functions 
export add_axes!, add_axes_projected!, add_axes_rotating!, add_axes_fixedoffset!,
       add_point!, add_point_dynamical!, add_point_fixedoffset!,
       add_direction!, add_axes_alias!, add_point_alias!,
       add_point_ephemeris!, add_axes_ephemeris!

include("core/axes.jl")
include("core/points.jl")
include("core/directions.jl")

# Transformations
export rotation3, rotation6, rotation9, rotation12,
       vector3, vector6, vector9, vector12,
       direction3, direction6, direction9, direction12

include("core/transform.jl")

# Prepared storage path and fully specialized compiled fast-path
export PreparedRotation, PreparedTranslation, PreparedDirection,
       CompiledRotation, CompiledTranslation, CompiledDirection,
       prepare_rotation, prepare_translation, prepare_direction,
       compile_rotation, compile_translation, compile_direction

include("core/compiled.jl")

# ==========================================================================================
# Definitions
# ==========================================================================================

export AXESID_ICRF, AXESID_GCRF,
       AXESID_ECL2000, AXESID_EME2000, AXESID_TEME,
       AXESID_MOONME_DE421, AXESID_MOONPA_DE421, AXESID_MOONPA_DE440

include("definitions/index.jl")

export add_axes_icrf!, add_axes_gcrf!, add_axes_eme2000!, add_axes_ecl2000!

include("definitions/celestial.jl")
include("definitions/ecliptic.jl")

export add_point_ephemeris!

include("definitions/ephemeris.jl")

export add_axes_frozen!

include("definitions/frozen.jl")

export add_axes_itrf!, add_axes_cirf!, add_axes_tirf!,
       add_axes_mod!, add_axes_tod!, add_axes_gtod!, add_axes_pef!

include("definitions/terrestrial.jl")

export teme_rot3_gcrf_to_teme, teme_rot6_gcrf_to_teme,
       teme_rot9_gcrf_to_teme, teme_rot12_gcrf_to_teme, add_axes_teme!

include("definitions/teme.jl")

export add_axes_bci2000!, add_axes_bcrtod!

include("definitions/planetary.jl")

export add_axes_pa440!, add_axes_pa421!, add_axes_me421!

include("definitions/lunar.jl")

export add_axes_topocentric!, add_point_surface!

include("definitions/topocentric.jl")

export add_direction_position!, add_direction_velocity!, add_direction_orthogonal!,
       add_direction_fixed!

include("definitions/directions.jl")

export add_axes_twodir!

include("definitions/axesfromdir.jl")

export add_axes_fixed_quaternion!, add_axes_fixed_angles!, add_axes_fixed_angleaxis!

include("definitions/attitude.jl")

end
