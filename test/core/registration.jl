using FrameTransformations
using ForwardDiff
using ReferenceFrameRotations
using StaticArrays
using Test

@testset "Cumulative derivative registration" begin
    @testset "position-only operations do not touch higher providers" begin
        frames = FrameSystem{4,Float64}()
        add_axes!(frames, :I, 1)
        add_axes_rotating!(
            frames, :A, 2, 1, t -> angle_to_dcm(t, :Z);
            rotation6=t -> error("rotation6 must not be called"),
            rotation9=t -> error("rotation9 must not be called"),
            rotation12=t -> error("rotation12 must not be called"),
        )
        add_point!(frames, :O, 1, 1)
        add_point_dynamical!(
            frames, :P, 2, 1, 1, t -> SVector(t, t^2, one(t));
            state6=t -> error("state6 must not be called"),
            state9=t -> error("state9 must not be called"),
            state12=t -> error("state12 must not be called"),
        )
        add_direction!(
            frames, :D, 1, t -> SVector(cos(t), sin(t), zero(t));
            direction6=t -> error("direction6 must not be called"),
            direction9=t -> error("direction9 must not be called"),
            direction12=t -> error("direction12 must not be called"),
        )

        for rotation in (
            t -> rotation3(frames, :I, :A, t),
            prepare_rotation(frames, :I, :A, Val(1)),
            compile_rotation(frames, :I, :A, Val(1)),
        )
            @test rotation(0.2) isa Rotation{1}
        end
        for translation in (
            t -> vector3(frames, :O, :P, :I, t),
            prepare_translation(frames, :O, :P, :I, Val(1)),
            compile_translation(frames, :O, :P, :I, Val(1)),
        )
            @test translation(0.2) isa SVector{3}
        end
        for direction in (
            t -> direction3(frames, :D, :I, t),
            prepare_direction(frames, :D, :I, Val(1)),
            compile_direction(frames, :D, :I, Val(1)),
        )
            @test direction(0.2) isa SVector{3}
        end
    end

    @testset "missing orders derive from the highest cumulative provider" begin
        frames = FrameSystem{3,Float64}()
        add_axes!(frames, :I, 1)

        rotation3_function(t) = angle_to_dcm(t^2, :Z)
        rotation6_calls = Ref(0)
        rotation6_function(t) = begin
            rotation6_calls[] += 1
            Rotation(
                rotation3_function(t),
                ForwardDiff.derivative(rotation3_function, t),
            )
        end
        add_axes_rotating!(
            frames, :A, 2, 1, rotation3_function;
            rotation6=rotation6_function,
        )

        position_calls = Ref(0)
        position(t) = begin
            position_calls[] += 1
            SVector(t, t^2, t^3)
        end
        state6_calls = Ref(0)
        state6(t) = begin
            state6_calls[] += 1
            SVector(t, t^2, t^3, one(t), 2t, 3t^2)
        end
        add_point!(frames, :O, 1, 1)
        add_point_dynamical!(frames, :P, 2, 1, 1, position; state6)

        direction3_calls = Ref(0)
        direction3_function(t) = begin
            direction3_calls[] += 1
            SVector(t, t^2, t^3)
        end
        direction6_calls = Ref(0)
        direction6_function(t) = begin
            direction6_calls[] += 1
            SVector(t, t^2, t^3, one(t), 2t, 3t^2)
        end
        add_direction!(
            frames, :D, 1, direction3_function;
            direction6=direction6_function,
        )

        time = 0.4
        rotation = rotation9(frames, :I, :A, time)
        state = vector9(frames, :O, :P, :I, time)
        direction = direction9(frames, :D, :I, time)

        @test rotation6_calls[] > 0
        @test rotation[3] ≈ ForwardDiff.derivative(
            epoch -> rotation6_function(epoch)[2], time)
        @test position_calls[] == 0
        @test state6_calls[] > 0
        @test state[7:9] ≈ SVector(0.0, 2.0, 6time)
        @test direction3_calls[] == 0
        @test direction6_calls[] > 0
        @test direction[7:9] ≈ SVector(0.0, 2.0, 6time)
    end

    @testset "explicit cumulative orders take priority" begin
        frames = FrameSystem{3,Float64}()
        add_axes!(frames, :I, 1)
        explicit_rotation9(t) = Rotation(
            angle_to_dcm(t, :Z),
            DCM(zero(t) * I),
            DCM(one(t) * I),
        )
        add_axes_rotating!(
            frames, :A, 2, 1, t -> angle_to_dcm(t, :Z);
            rotation6=t -> error("rotation6 must not be called"),
            rotation9=explicit_rotation9,
        )

        add_point!(frames, :O, 1, 1)
        explicit_state9(t) = SVector(
            t, t, t, zero(t), zero(t), zero(t), one(t), one(t), one(t))
        add_point_dynamical!(
            frames, :P, 2, 1, 1, t -> error("position must not be called");
            state6=t -> error("state6 must not be called"),
            state9=explicit_state9,
        )
        add_direction!(
            frames, :D, 1, t -> error("direction3 must not be called");
            direction6=t -> error("direction6 must not be called"),
            direction9=explicit_state9,
        )

        @test rotation9(frames, :I, :A, 0.2)[3] ≈ DCM(1.0I)
        @test vector9(frames, :O, :P, :I, 0.2) == explicit_state9(0.2)
        @test direction9(frames, :D, :I, 0.2) == explicit_state9(0.2)
    end
end
