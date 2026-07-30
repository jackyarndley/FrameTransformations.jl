using FrameTransformations
using DifferentiationInterface
using ReferenceFrameRotations
using StaticArrays
using ForwardDiff
using FiniteDiff
using ChainRulesCore
using Zygote
using Mooncake
using Test

# Build a FrameSystem with analytical point and rotating axes for testing
function build_test_framesystem()
    fr = FrameSystem{4,Float64}()
    add_axes!(fr, :ICRF, 1)

    # Simple rotating axes around Z
    ω = 0.001
    add_axes_rotating!(fr, :RotZ, 2, 1, t -> angle_to_dcm(ω * t, :Z))

    # Root point + dynamical point with known analytical position/velocity
    add_point!(fr, :Root, 1, 1)
    add_point_dynamical!(fr, :Target, 2, 1, 1, t -> SVector(cos(t), sin(t), 1.0))
    add_direction!(fr, :LineOfSight, 1, t -> SVector(cos(t), sin(t), 0.0))

    return fr
end

const _fr = build_test_framesystem()

@testset "DifferentiationInterface backend compatibility" begin
    time = 1.0
    forward_backend = AutoForwardDiff()
    finite_backend = AutoFiniteDiff(fdtype=Val(:central))
    reverse_backend = AutoZygote()

    position = epoch -> vector3(_fr, 1, 2, 1, epoch)
    state6 = epoch -> vector6(_fr, 1, 2, 1, epoch)
    rotation_sum = epoch -> sum(rotation3(_fr, 1, 2, epoch).m[1])

    expected_velocity = vector6(_fr, 1, 2, 1, time)[4:6]
    expected_state6_derivative = vector9(_fr, 1, 2, 1, time)[4:9]
    expected_rotation_sum_derivative =
        sum(rotation6(_fr, 1, 2, time).m[2])

    for backend in (forward_backend, reverse_backend)
        @test derivative(position, backend, time) ≈
              expected_velocity atol=1e-12 rtol=1e-12
        @test derivative(state6, backend, time) ≈
              expected_state6_derivative atol=1e-12 rtol=1e-12
        @test derivative(rotation_sum, backend, time) ≈
              expected_rotation_sum_derivative atol=1e-12 rtol=1e-12
    end

    @test derivative(position, finite_backend, time) ≈
          expected_velocity atol=1e-7 rtol=1e-7
    @test derivative(state6, finite_backend, time) ≈
          expected_state6_derivative atol=1e-6 rtol=1e-6
end

@testset "analytic rule coverage for every transform order" begin
    time = 0.7
    reverse_backend = AutoZygote()

    for (function_object, next_function) in (
        (vector3, vector6),
        (vector6, vector9),
        (vector9, vector12),
    )
        scalar_state = epoch -> sum(function_object(_fr, 1, 2, 1, epoch))
        next_state = next_function(_fr, 1, 2, 1, time)
        width = length(function_object(_fr, 1, 2, 1, time))
        expected = sum(next_state[4:(3 + width)])
        @test derivative(scalar_state, reverse_backend, time) ≈
              expected atol=1e-12 rtol=1e-12
    end

    for (function_object, next_function) in (
        (direction3, direction6),
        (direction6, direction9),
        (direction9, direction12),
    )
        scalar_direction =
            epoch -> sum(function_object(_fr, :LineOfSight, 1, epoch))
        next_state = next_function(_fr, :LineOfSight, 1, time)
        width = length(function_object(_fr, :LineOfSight, 1, time))
        expected = sum(next_state[4:(3 + width)])
        @test derivative(scalar_direction, reverse_backend, time) ≈
              expected atol=1e-12 rtol=1e-12
    end

    for (function_object, next_function) in (
        (rotation3, rotation6),
        (rotation6, rotation9),
        (rotation9, rotation12),
    )
        scalar_rotation =
            epoch -> sum(sum, function_object(_fr, 1, 2, epoch).m)
        next_state = next_function(_fr, 1, 2, time)
        component_count = length(function_object(_fr, 1, 2, time))
        expected = sum(sum, next_state.m[2:(component_count + 1)])
        @test derivative(scalar_rotation, reverse_backend, time) ≈
              expected atol=1e-12 rtol=1e-12
    end

    for scalar_function in (
        epoch -> sum(vector12(_fr, 1, 2, 1, epoch)),
        epoch -> sum(direction12(_fr, :LineOfSight, 1, epoch)),
        epoch -> sum(sum, rotation12(_fr, 1, 2, epoch).m),
    )
        expected = ForwardDiff.derivative(scalar_function, time)
        @test derivative(scalar_function, reverse_backend, time) ≈
              expected atol=1e-11 rtol=1e-11
    end
end

# --------------------------------------------------------------------------
# vector3
# --------------------------------------------------------------------------
@testset "vector3 — Zygote vs ForwardDiff" begin
    for t0 in (0.0, 1.0, π / 3)
        fd = ForwardDiff.derivative(t -> sum(vector3(_fr, 1, 2, 1, t)), t0)
        zy = Zygote.gradient(t -> sum(vector3(_fr, 1, 2, 1, t)), t0)[1]
        @test fd ≈ zy rtol = 1e-10
    end
end

@testset "vector3 — Mooncake vs ForwardDiff" begin
    for t0 in (0.0, 1.0, π / 3)
        fd = ForwardDiff.derivative(t -> sum(vector3(_fr, 1, 2, 1, t)), t0)
        f = t -> sum(vector3(_fr, 1, 2, 1, t))
        rule = Mooncake.build_rrule(f, t0)
        _, (_, mc) = Mooncake.value_and_gradient!!(rule, f, t0)
        @test fd ≈ mc rtol = 1e-10
    end
end

# --------------------------------------------------------------------------
# vector6
# --------------------------------------------------------------------------
@testset "vector6 — Zygote vs ForwardDiff" begin
    for t0 in (0.0, 1.0, π / 3)
        fd = ForwardDiff.derivative(t -> sum(vector6(_fr, 1, 2, 1, t)), t0)
        zy = Zygote.gradient(t -> sum(vector6(_fr, 1, 2, 1, t)), t0)[1]
        @test fd ≈ zy rtol = 1e-10
    end
end

@testset "vector6 — Mooncake vs ForwardDiff" begin
    for t0 in (0.0, 1.0, π / 3)
        fd = ForwardDiff.derivative(t -> sum(vector6(_fr, 1, 2, 1, t)), t0)
        f = t -> sum(vector6(_fr, 1, 2, 1, t))
        rule = Mooncake.build_rrule(f, t0)
        _, (_, mc) = Mooncake.value_and_gradient!!(rule, f, t0)
        @test fd ≈ mc rtol = 1e-10
    end
end

# --------------------------------------------------------------------------
# rotation3 — access .m[1] (DCM) which is how downstream code uses it
# --------------------------------------------------------------------------
@testset "rotation3 — Zygote vs ForwardDiff" begin
    for t0 in (0.0, 1.0, π / 3)
        fd = ForwardDiff.derivative(t -> sum(rotation3(_fr, 1, 2, t).m[1]), t0)
        zy = Zygote.gradient(t -> sum(rotation3(_fr, 1, 2, t).m[1]), t0)[1]
        @test fd ≈ zy rtol = 1e-10
    end
end

# --------------------------------------------------------------------------
# rrule correctness: pullback matches FiniteDiff
# --------------------------------------------------------------------------
@testset "rrule pullback vs FiniteDiff" begin
    t0 = 1.0

    # vector3
    g_fd = FiniteDiff.finite_difference_derivative(
        t -> sum(vector3(_fr, 1, 2, 1, t)), t0
    )
    g_zy = Zygote.gradient(t -> sum(vector3(_fr, 1, 2, 1, t)), t0)[1]
    @test g_fd ≈ g_zy rtol = 1e-5

    # vector6
    g_fd6 = FiniteDiff.finite_difference_derivative(
        t -> sum(vector6(_fr, 1, 2, 1, t)), t0
    )
    g_zy6 = Zygote.gradient(t -> sum(vector6(_fr, 1, 2, 1, t)), t0)[1]
    @test g_fd6 ≈ g_zy6 rtol = 1e-5

    # rotation3
    g_fd_r = FiniteDiff.finite_difference_derivative(
        t -> sum(rotation3(_fr, 1, 2, t).m[1]), t0
    )
    g_zy_r = Zygote.gradient(t -> sum(rotation3(_fr, 1, 2, t).m[1]), t0)[1]
    @test g_fd_r ≈ g_zy_r rtol = 1e-5
end
