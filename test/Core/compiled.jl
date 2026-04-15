using FrameTransformations
using ReferenceFrameRotations
using StaticArrays
using ForwardDiff
using Test

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
    # Rotation: 2-node (direct parent-child)
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
    # Rotation: inverse direction
    # --------------------------------------------------------------------------
    @testset "compile_rotation — inverse" begin
        cr_fwd = compile_rotation(fr, :ICRF, :A)
        cr_inv = compile_rotation(fr, :A, :ICRF)
        @test cr_fwd(t)[1] ≈ inv(cr_inv(t))[1]
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
    # ForwardDiff through compiled rotation
    # --------------------------------------------------------------------------
    @testset "ForwardDiff — compiled rotation" begin
        cr = compile_rotation(fr, :ICRF, :A)
        dcm_at_t(t) = cr(t)[1]
        J = ForwardDiff.derivative(dcm_at_t, t)
        R = cr(t)
        @test J ≈ R[2] atol=1e-10
    end

    # --------------------------------------------------------------------------
    # ForwardDiff through compiled translation
    # --------------------------------------------------------------------------
    @testset "ForwardDiff — compiled translation" begin
        ct = compile_translation(fr, :Origin, :P1, :ICRF)
        pos(t) = ct(t)[SVector(1,2,3)]
        J = ForwardDiff.derivative(pos, t)
        v = ct(t)
        @test J ≈ v[SVector(4,5,6)] atol=1e-10
    end

    # --------------------------------------------------------------------------
    # ForwardDiff through compiled direction
    # --------------------------------------------------------------------------
    @testset "ForwardDiff — compiled direction" begin
        cd = compile_direction(fr, :sun, :ICRF)
        dir_pos(t) = cd(t)[SVector(1,2,3)]
        J = ForwardDiff.derivative(dir_pos, t)
        d = cd(t)
        @test J ≈ d[SVector(4,5,6)] atol=1e-10
    end

end
