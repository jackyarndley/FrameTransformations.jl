using FrameTransformations
using Ephemerides: EphemerisProvider
using ReferenceFrameRotations: angle_to_dcm
using StaticArrays: SVector

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const TEST_ASSETS = joinpath(REPO_ROOT, "test", "assets")

asset_path(name) = joinpath(TEST_ASSETS, name)

function rotation_direct_callable(fr, from, to, order)
    if order == 2
        return t -> rotation6(fr, from, to, t)
    elseif order == 3
        return t -> rotation9(fr, from, to, t)
    elseif order == 4
        return t -> rotation12(fr, from, to, t)
    end
    throw(ArgumentError("unsupported rotation order $order"))
end

function rotation_compiled_callable(fr, from, to, order)
    if order == 2
        return compile_rotation(fr, from, to, Val(2))
    elseif order == 3
        return compile_rotation(fr, from, to, Val(3))
    elseif order == 4
        return compile_rotation(fr, from, to, Val(4))
    end
    throw(ArgumentError("unsupported rotation order $order"))
end

function rotation_prepared_callable(fr, from, to, order)
    return prepare_rotation(fr, from, to, Val(order))
end

function translation_direct_callable(fr, from, to, axes, order)
    if order == 2
        return t -> vector6(fr, from, to, axes, t)
    elseif order == 3
        return t -> vector9(fr, from, to, axes, t)
    elseif order == 4
        return t -> vector12(fr, from, to, axes, t)
    end
    throw(ArgumentError("unsupported translation order $order"))
end

function translation_compiled_callable(fr, from, to, axes, order)
    if order == 2
        return compile_translation(fr, from, to, axes, Val(2))
    elseif order == 3
        return compile_translation(fr, from, to, axes, Val(3))
    elseif order == 4
        return compile_translation(fr, from, to, axes, Val(4))
    end
    throw(ArgumentError("unsupported translation order $order"))
end

function translation_prepared_callable(fr, from, to, axes, order)
    return prepare_translation(fr, from, to, axes, Val(order))
end

function direction_direct_callable(fr, name, axes, order)
    if order == 2
        return t -> direction6(fr, name, axes, t)
    elseif order == 3
        return t -> direction9(fr, name, axes, t)
    elseif order == 4
        return t -> direction12(fr, name, axes, t)
    end
    throw(ArgumentError("unsupported direction order $order"))
end

function direction_compiled_callable(fr, name, axes, order)
    if order == 2
        return compile_direction(fr, name, axes, Val(2))
    elseif order == 3
        return compile_direction(fr, name, axes, Val(3))
    elseif order == 4
        return compile_direction(fr, name, axes, Val(4))
    end
    throw(ArgumentError("unsupported direction order $order"))
end

function direction_prepared_callable(fr, name, axes, order)
    return prepare_direction(fr, name, axes, Val(order))
end

function order_suffix(order)
    return 3 * order
end

function build_synthetic_benchmark_scenario(order::Int = 4)
    fr = FrameSystem{order,Float64}()
    add_axes!(fr, :ICRF, 1)

    add_axes_rotating!(fr, :A, 2, 1, t -> angle_to_dcm(t, :Z))
    add_axes_rotating!(fr, :B, 3, 2, t -> angle_to_dcm(0.5 * t, :X))
    add_axes_rotating!(fr, :C, 4, 3, t -> angle_to_dcm(0.3 * t, :Y))

    add_point!(fr, :Origin, 1, 1)
    add_point_dynamical!(fr, :P1, 2, 1, 1, t -> SVector(cos(t), sin(t), 1.0))
    add_point_dynamical!(fr, :P2, 3, 2, 2, t -> SVector(t^2, 2 * t, 3.0))

    add_direction!(fr, :sun, 1, t -> SVector(cos(t), sin(t), 0.0))

    t0 = pi / 3
    rotation_direct = rotation_direct_callable(fr, :ICRF, :C, order)
    rotation_prepared = rotation_prepared_callable(fr, :ICRF, :C, order)
    rotation_compiled = rotation_compiled_callable(fr, :ICRF, :C, order)
    translation_direct = translation_direct_callable(fr, :Origin, :P2, :A, order)
    translation_prepared = translation_prepared_callable(
        fr, :Origin, :P2, :A, order)
    translation_compiled = translation_compiled_callable(fr, :Origin, :P2, :A, order)
    direction_direct = direction_direct_callable(fr, :sun, :A, order)
    direction_prepared = direction_prepared_callable(fr, :sun, :A, order)
    direction_compiled = direction_compiled_callable(fr, :sun, :A, order)

    return (
        order = order,
        suffix = order_suffix(order),
        fr = fr,
        t0 = t0,
        rotation_direct = rotation_direct,
        rotation_prepared = rotation_prepared,
        rotation_compiled = rotation_compiled,
        translation_direct = translation_direct,
        translation_prepared = translation_prepared,
        translation_compiled = translation_compiled,
        direction_direct = direction_direct,
        direction_prepared = direction_prepared,
        direction_compiled = direction_compiled,
    )
end

function build_de440_benchmark_scenario(order::Int = 4)
    fr = FrameSystem{order,Float64}()
    add_axes_icrf!(fr)

    eph = EphemerisProvider(asset_path("moon_pa_de440_200625.bpc"))
    add_axes_pa440!(fr, eph, :PA440)
    add_axes_me421!(fr, :ME421, :PA440)

    t0 = 1.0e6
    rotation_direct = rotation_direct_callable(fr, :ME421, :ICRF, order)
    rotation_prepared = rotation_prepared_callable(fr, :ME421, :ICRF, order)
    rotation_compiled = rotation_compiled_callable(fr, :ME421, :ICRF, order)

    @assert rotation_direct(t0) == rotation_compiled(t0)

    return (
        order = order,
        suffix = order_suffix(order),
        fr = fr,
        t0 = t0,
        rotation_direct = rotation_direct,
        rotation_prepared = rotation_prepared,
        rotation_compiled = rotation_compiled,
    )
end

function build_de440_vector_benchmark_scenario(order::Int = 4)
    fr = FrameSystem{order,Float64}()
    add_axes_icrf!(fr)

    lunar_eph = EphemerisProvider(asset_path("moon_pa_de440_200625.bpc"))
    point_eph = EphemerisProvider(asset_path("de432s.bsp"))

    add_axes_pa440!(fr, lunar_eph, :PA440)
    add_axes_me421!(fr, :ME421, :PA440)

    add_point!(fr, :SSB, 0, 1)
    add_point_ephemeris!(fr, point_eph, :EarthB, 3)
    add_point_ephemeris!(fr, point_eph, :Earth, 399)
    add_point_ephemeris!(fr, point_eph, :Moon, 301)

    t0 = 1.0e6
    translation_direct = translation_direct_callable(fr, :Earth, :Moon, :ME421, order)
    translation_prepared = translation_prepared_callable(
        fr, :Earth, :Moon, :ME421, order)
    translation_compiled = translation_compiled_callable(fr, :Earth, :Moon, :ME421, order)

    @assert translation_direct(t0) ≈ translation_compiled(t0)

    return (
        order = order,
        suffix = order_suffix(order),
        fr = fr,
        t0 = t0,
        translation_direct = translation_direct,
        translation_prepared = translation_prepared,
        translation_compiled = translation_compiled,
    )
end
