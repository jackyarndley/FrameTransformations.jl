using FrameTransformations
using DifferentiationInterface: AutoForwardDiff, derivative
using ReferenceFrameRotations
using StaticArrays
import ForwardDiff
using Test

const backend = AutoForwardDiff()

@testset "Compiled fast-path" begin

    # --------------------------------------------------------------------------
    # Setup: build a frame system with a 3-hop axes chain and multi-hop points
    # --------------------------------------------------------------------------
    fr = FrameSystem{4, Float64}()
    add_axes!(fr, :ICRF, 1)

    rA(t) = angle_to_dcm(t, :Z)
    add_axes_rotating!(fr, :A, 2, 1, rA)

    rB(t) = angle_to_dcm(0.5t, :X)
    add_axes_rotating!(fr, :B, 3, 2, rB)

    rC(t) = angle_to_dcm(0.3t, :Y)
    add_axes_rotating!(fr, :C, 4, 3, rC)

    add_point!(fr, :Origin, 1, 1)

    fP(t) = SVector(cos(t), sin(t), 1.0)
    add_point_dynamical!(fr, :P1, 2, 1, 1, fP)

    fQ(t) = SVector(t^2, 2t, 3.0)
    add_point_dynamical!(fr, :P2, 3, 2, 2, fQ)

    t = π / 3

    # --------------------------------------------------------------------------
    # Rotation: 2-node (direct parent-child), default order (O=4)
    # --------------------------------------------------------------------------
    @testset "compile_rotation — 2-node" begin
        cr = compile_rotation(fr, :ICRF, :A)
        @test cr isa CompiledRotation

        R_core = rotation12(fr, 1, 2, t)
        R_comp = cr(t)
        for i in 1:4
            @test R_comp[i] ≈ R_core[i]
        end
    end

    # --------------------------------------------------------------------------
    # Rotation: identity (same axes)
    # --------------------------------------------------------------------------
    @testset "compile_rotation — identity" begin
        cr = compile_rotation(fr, :ICRF, :ICRF)
        R = cr(t)
        @test R[1] ≈ one(DCM{Float64})
    end

    # --------------------------------------------------------------------------
    # Rotation: inverse direction is part of the specialized operation
    # --------------------------------------------------------------------------
    @testset "compile_rotation — inverse" begin
        cr_fwd = compile_rotation(fr, :ICRF, :A)
        cr_inv = compile_rotation(fr, :A, :ICRF)
        @test cr_fwd(t)[1] ≈ inv(cr_inv(t))[1]

        @test cr_fwd isa CompiledRotation{4}
        @test cr_inv isa CompiledRotation{4}
        @test typeof(cr_fwd) !== typeof(cr_inv)
    end

    # --------------------------------------------------------------------------
    # Rotation: multi-hop (ICRF → A → B → C)
    # --------------------------------------------------------------------------
    @testset "compile_rotation — multi-hop" begin
        cr = compile_rotation(fr, :ICRF, :C)
        R_core = rotation12(fr, 1, 4, t)
        R_comp = cr(t)
        for i in 1:4
            @test R_comp[i] ≈ R_core[i]
        end

        cr_ab = compile_rotation(fr, :A, :B)
        R_ab_core = rotation12(fr, 2, 3, t)
        R_ab_comp = cr_ab(t)
        for i in 1:4
            @test R_ab_comp[i] ≈ R_ab_core[i]
        end
    end

    # --------------------------------------------------------------------------
    # Rotation: reverse multi-hop (C → B → A → ICRF)
    # --------------------------------------------------------------------------
    @testset "compile_rotation — reverse multi-hop" begin
        cr_rev = compile_rotation(fr, :C, :ICRF)
        R_core = rotation12(fr, 4, 1, t)
        R_comp = cr_rev(t)
        for i in 1:4
            @test R_comp[i] ≈ R_core[i]
        end

        cr_fwd = compile_rotation(fr, :ICRF, :C)
        @test cr_rev(t)[1] ≈ inv(cr_fwd(t))[1]
    end

    # --------------------------------------------------------------------------
    # Rotation: order selection — Val{1} gives only DCM, no derivatives
    # --------------------------------------------------------------------------
    @testset "compile_rotation — order selection" begin
        cr1 = compile_rotation(fr, :ICRF, :A, Val(1))
        @test cr1 isa CompiledRotation{1}
        @test cr1.function_object(t) isa Rotation{1}
        R1 = cr1(t)
        @test R1[1] ≈ rotation3(fr, 1, 2, t)[1]

        cr2 = compile_rotation(fr, :ICRF, :A, Val(2))
        @test cr2 isa CompiledRotation{2}
        @test cr2.function_object(t) isa Rotation{2}
        R2 = cr2(t)
        R2_core = rotation6(fr, 1, 2, t)
        @test R2[1] ≈ R2_core[1]
        @test R2[2] ≈ R2_core[2]

        # Order exceeding system order should error
        @test_throws ArgumentError compile_rotation(fr, :ICRF, :A, Val(5))

        # Multi-hop with reduced order
        cr1_mh = compile_rotation(fr, :ICRF, :C, Val(1))
        R1_mh = cr1_mh(t)
        @test R1_mh[1] ≈ rotation3(fr, 1, 4, t)[1]
    end

    # --------------------------------------------------------------------------
    # Translation: 2-node, same axes
    # --------------------------------------------------------------------------
    @testset "compile_translation — 2-node same axes" begin
        ct = compile_translation(fr, :Origin, :P1, :ICRF)
        @test ct isa CompiledTranslation

        v_core = vector12(fr, 1, 2, 1, t)
        v_comp = ct(t)
        @test v_comp ≈ v_core
    end

    # --------------------------------------------------------------------------
    # Translation: 2-node, different axes (requires rotation)
    # --------------------------------------------------------------------------
    @testset "compile_translation — 2-node different axes" begin
        ct = compile_translation(fr, :Origin, :P1, :A)
        v_core = vector12(fr, 1, 2, 2, t)
        v_comp = ct(t)
        @test v_comp ≈ v_core
    end

    # --------------------------------------------------------------------------
    # Translation: identity (same point)
    # --------------------------------------------------------------------------
    @testset "compile_translation — identity" begin
        ct = compile_translation(fr, :Origin, :Origin, :ICRF)
        @test ct(t) ≈ @SVector zeros(12)
    end

    # --------------------------------------------------------------------------
    # Translation: multi-hop (Origin → P1 → P2)
    # --------------------------------------------------------------------------
    @testset "compile_translation — multi-hop" begin
        ct = compile_translation(fr, :Origin, :P2, :ICRF)
        v_core = vector12(fr, 1, 3, 1, t)
        v_comp = ct(t)
        @test v_comp ≈ v_core
    end

    # --------------------------------------------------------------------------
    # Translation: multi-hop in a different axes
    # --------------------------------------------------------------------------
    @testset "compile_translation — multi-hop different axes" begin
        ct = compile_translation(fr, :Origin, :P2, :A)
        v_core = vector12(fr, 1, 3, 2, t)
        v_comp = ct(t)
        @test v_comp ≈ v_core
    end

    # --------------------------------------------------------------------------
    # Translation: reverse multi-hop (P2 → P1 → Origin), exercises inv_flag
    # in _compile_translation_forward's loop body
    # --------------------------------------------------------------------------
    @testset "compile_translation — reverse multi-hop" begin
        ct_rev = compile_translation(fr, :P2, :Origin, :ICRF)
        v_core = vector12(fr, 3, 1, 1, t)
        @test ct_rev(t) ≈ v_core

        ct_rev_a = compile_translation(fr, :P2, :Origin, :A)
        v_core_a = vector12(fr, 3, 1, 2, t)
        @test ct_rev_a(t) ≈ v_core_a
    end

    # --------------------------------------------------------------------------
    # Translation: order selection
    # --------------------------------------------------------------------------
    @testset "compile_translation — order selection" begin
        ct1 = compile_translation(fr, :Origin, :P1, :ICRF, Val(1))
        @test ct1 isa CompiledTranslation{1}
        v1 = ct1(t)
        @test v1 ≈ vector3(fr, 1, 2, 1, t)

        ct2 = compile_translation(fr, :Origin, :P1, :A, Val(2))
        @test ct2 isa CompiledTranslation{2}
        v2 = ct2(t)
        @test v2 ≈ vector6(fr, 1, 2, 2, t)

        @test_throws ArgumentError compile_translation(fr, :Origin, :P1, :ICRF, Val(5))
    end

    # --------------------------------------------------------------------------
    # Direction
    # --------------------------------------------------------------------------
    @testset "compile_direction" begin
        dfun(t) = SVector(cos(t), sin(t), 0.0)
        add_direction!(fr, :sun, 1, dfun)

        cd = compile_direction(fr, :sun, :ICRF)
        @test cd isa CompiledDirection
        d_core = direction12(fr, :sun, 1, t)
        d_comp = cd(t)
        @test d_comp ≈ d_core

        # Direction in different axes
        cd2 = compile_direction(fr, :sun, :A)
        d_core2 = direction12(fr, :sun, 2, t)
        d_comp2 = cd2(t)
        @test d_comp2 ≈ d_core2
    end

    # --------------------------------------------------------------------------
    # Direction: order selection
    # --------------------------------------------------------------------------
    @testset "compile_direction — order selection" begin
        cd1 = compile_direction(fr, :sun, :ICRF, Val(1))
        @test cd1 isa CompiledDirection{1}
        d1 = cd1(t)
        @test d1 ≈ direction3(fr, :sun, 1, t)

        @test_throws ArgumentError compile_direction(fr, :sun, :ICRF, Val(5))
    end

    @testset "direct, prepared, and compiled equivalence" begin
        rotation_functions = (rotation3, rotation6, rotation9, rotation12)
        vector_functions = (vector3, vector6, vector9, vector12)
        direction_functions = (direction3, direction6, direction9, direction12)

        for index in 1:4
            order_value = Val(index)

            for (from, to) in (
                (:ICRF, :ICRF), # identity
                (:ICRF, :A),    # direct forward
                (:A, :ICRF),    # direct inverse
                (:ICRF, :C),    # multi-hop forward
                (:C, :ICRF),    # multi-hop reverse
            )
                direct = rotation_functions[index](fr, from, to, t)
                prepared = prepare_rotation(fr, from, to, order_value)(t)
                compiled = compile_rotation(fr, from, to, order_value)(t)
                for component in 1:index
                    @test prepared[component] ≈ direct[component]
                    @test compiled[component] ≈ direct[component]
                end
            end

            for (from, to, axes) in (
                (:Origin, :Origin, :ICRF),
                (:Origin, :P1, :ICRF),
                (:P1, :Origin, :A),
                (:Origin, :P2, :C),
                (:P2, :Origin, :A),
            )
                direct = vector_functions[index](fr, from, to, axes, t)
                prepared = prepare_translation(
                    fr, from, to, axes, order_value)(t)
                compiled = compile_translation(
                    fr, from, to, axes, order_value)(t)
                @test prepared ≈ direct
                @test compiled ≈ direct
                @test prepared isa SVector{3 * index}
            end

            for axes in (:ICRF, :A, :C)
                direct = direction_functions[index](fr, :sun, axes, t)
                prepared = prepare_direction(fr, :sun, axes, order_value)(t)
                compiled = compile_direction(fr, :sun, axes, order_value)(t)
                @test prepared ≈ direct
                @test compiled ≈ direct
                @test prepared isa SVector{3 * index}
            end
        end
    end

    # --------------------------------------------------------------------------
    # DifferentiationInterface through compiled rotation
    # --------------------------------------------------------------------------
    @testset "DifferentiationInterface — compiled rotation" begin
        cr = compile_rotation(fr, :ICRF, :A)
        dcm_at_t(t) = cr(t)[1]
        J = derivative(dcm_at_t, backend, t)
        R = cr(t)
        @test J ≈ R[2] atol=1e-10
    end

    # --------------------------------------------------------------------------
    # DifferentiationInterface through compiled translation
    # --------------------------------------------------------------------------
    @testset "DifferentiationInterface — compiled translation" begin
        ct = compile_translation(fr, :Origin, :P1, :ICRF)
        pos(t) = ct(t)[SVector(1,2,3)]
        J = derivative(pos, backend, t)
        v = ct(t)
        @test J ≈ v[SVector(4,5,6)] atol=1e-10
    end

    # --------------------------------------------------------------------------
    # DifferentiationInterface through compiled direction
    # --------------------------------------------------------------------------
    @testset "DifferentiationInterface — compiled direction" begin
        cd = compile_direction(fr, :sun, :ICRF)
        dir_pos(t) = cd(t)[SVector(1,2,3)]
        J = derivative(dir_pos, backend, t)
        d = cd(t)
        @test J ≈ d[SVector(4,5,6)] atol=1e-10
    end

    @testset "DifferentiationInterface — prepared routes" begin
        prepared_rotation = prepare_rotation(fr, :C, :ICRF, Val(2))
        prepared_translation = prepare_translation(
            fr, :P2, :Origin, :A, Val(2))
        prepared_direction = prepare_direction(fr, :sun, :C, Val(2))

        rotation_derivative = derivative(
            time -> prepared_rotation(time)[1], backend, t)
        translation_derivative = derivative(
            time -> prepared_translation(time)[SVector(1, 2, 3)], backend, t)
        direction_derivative = derivative(
            time -> prepared_direction(time)[SVector(1, 2, 3)], backend, t)

        @test rotation_derivative ≈ prepared_rotation(t)[2] atol=1e-10
        @test translation_derivative ≈
            prepared_translation(t)[SVector(4, 5, 6)] atol=1e-10
        @test direction_derivative ≈
            prepared_direction(t)[SVector(4, 5, 6)] atol=1e-10

        identity_rotation = prepare_rotation(fr, :ICRF, :ICRF, Val(1))
        identity_translation = prepare_translation(
            fr, :Origin, :Origin, :ICRF, Val(1))
        @test derivative(time -> identity_rotation(time)[1][1, 1], backend, t) == 0
        @test derivative(time -> identity_translation(time)[1], backend, t) == 0
    end

    # --------------------------------------------------------------------------
    # Prepared callables erase route structure from their public type
    # --------------------------------------------------------------------------
    @testset "prepared callables" begin
        compact_rotation_a = prepare_rotation(fr, :ICRF, :A, Val(1))
        compact_rotation_c = prepare_rotation(fr, :ICRF, :C, Val(1))
        compact_rotation_inverse = prepare_rotation(fr, :A, :ICRF, Val(1))
        compact_rotation_identity = prepare_rotation(fr, :ICRF, :ICRF, Val(1))
        @test typeof(compact_rotation_a) === typeof(compact_rotation_c)
        @test typeof(compact_rotation_a) === typeof(compact_rotation_inverse)
        @test typeof(compact_rotation_a) === typeof(compact_rotation_identity)
        @test compact_rotation_a isa PreparedRotation{1,Float64}
        @test compact_rotation_a(t)[1] ≈ rotation3(fr, 1, 2, t)[1]
        @test compact_rotation_c(t)[1] ≈ rotation3(fr, 1, 4, t)[1]
        @test compact_rotation_inverse(t)[1] ≈ rotation3(fr, 2, 1, t)[1]

        compact_translation_1 = prepare_translation(
            fr, :Origin, :P1, :ICRF, Val(1))
        compact_translation_2 = prepare_translation(
            fr, :Origin, :P2, :ICRF, Val(1))
        compact_translation_inverse = prepare_translation(
            fr, :P1, :Origin, :ICRF, Val(1))
        compact_translation_identity = prepare_translation(
            fr, :Origin, :Origin, :ICRF, Val(1))
        @test typeof(compact_translation_1) === typeof(compact_translation_2)
        @test typeof(compact_translation_1) === typeof(compact_translation_inverse)
        @test typeof(compact_translation_1) === typeof(compact_translation_identity)
        @test compact_translation_1 isa PreparedTranslation{1,Float64}
        @test compact_translation_1(t) ≈ vector3(fr, 1, 2, 1, t)
        @test compact_translation_2(t) ≈ vector3(fr, 1, 3, 1, t)

        compact_direction = prepare_direction(fr, :sun, :ICRF, Val(1))
        compact_direction_rotated = prepare_direction(fr, :sun, :C, Val(1))
        @test typeof(compact_direction) === typeof(compact_direction_rotated)
        @test compact_direction isa PreparedDirection{1,Float64}
        @test compact_direction(t) ≈ direction3(fr, :sun, 1, t)

        compact_translation_ad = prepare_translation(
            fr, :Origin, :P1, :ICRF, Val(2))
        compact_position(time) = compact_translation_ad(time)[SVector(1,2,3)]
        compact_jacobian = derivative(compact_position, backend, t)
        @test compact_jacobian ≈
            compact_translation_ad(t)[SVector(4,5,6)] atol=1e-10
    end

end
