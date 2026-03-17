using FrameTransformations
using ReferenceFrameRotations
using StaticArrays
using ForwardDiff
using Test

@testset "CompiledRotation" begin
    fr = FrameSystem{4, Float64}()
    add_axes!(fr, :ICRF, 1)

    r(t) = angle_to_dcm(t, :Z)
    add_axes_rotating!(fr, :Test, 2, 1, r)

    x = π / 3

    # --- Order 1 (rotation3) ---
    cr3 = compile_rotation3(fr, :ICRF, :Test)
    @test cr3(x).m[1] ≈ rotation3(fr, :ICRF, :Test, x).m[1]

    cr3_inv = compile_rotation3(fr, :Test, :ICRF)
    @test cr3_inv(x).m[1] ≈ rotation3(fr, :Test, :ICRF, x).m[1]

    # --- Order 2 (rotation6) ---
    cr6 = compile_rotation6(fr, :ICRF, :Test)
    R6c = cr6(x); R6r = rotation6(fr, :ICRF, :Test, x)
    @test R6c.m[1] ≈ R6r.m[1]
    @test R6c.m[2] ≈ R6r.m[2]

    cr6_inv = compile_rotation6(fr, :Test, :ICRF)
    R6ci = cr6_inv(x); R6ri = rotation6(fr, :Test, :ICRF, x)
    @test R6ci.m[1] ≈ R6ri.m[1]
    @test R6ci.m[2] ≈ R6ri.m[2]

    # --- Order 3 (rotation9) ---
    cr9 = compile_rotation9(fr, 1, 2)
    R9c = cr9(x); R9r = rotation9(fr, 1, 2, x)
    @test R9c.m[1] ≈ R9r.m[1]
    @test R9c.m[2] ≈ R9r.m[2]
    @test R9c.m[3] ≈ R9r.m[3]

    # --- Order 4 (rotation12) ---
    cr12 = compile_rotation12(fr, 1, 2)
    R12c = cr12(x); R12r = rotation12(fr, 1, 2, x)
    for i in 1:4
        @test R12c.m[i] ≈ R12r.m[i]
    end

    # Multiple time values
    for t_val in [0.0, 1.0, π/4, 2π]
        @test compile_rotation3(fr, 1, 2)(t_val).m[1] ≈ rotation3(fr, 1, 2, t_val).m[1]
    end
end

@testset "CompiledTranslation" begin
    fr = FrameSystem{4, Float64}()
    add_axes!(fr, :ICRF, 1)

    f(t) = SVector(cos(t), sin(t), 1.0)
    df(t) = ForwardDiff.derivative(f, t)

    add_point!(fr, :ROOT, 1, 1)
    add_point_dynamical!(fr, :Test, 2, 1, 1, f)

    x = π / 3

    # --- Order 1 (vector3) ---
    ct3 = compile_vector3(fr, 1, 2, 1)
    tr3 = ct3(x)
    v3r = vector3(fr, 1, 2, 1, x)
    @test SVector{3}(tr3[1][1], tr3[1][2], tr3[1][3]) ≈ v3r

    ct3_inv = compile_vector3(fr, 2, 1, 1)
    tr3_inv = ct3_inv(x)
    v3ri = vector3(fr, 2, 1, 1, x)
    @test SVector{3}(tr3_inv[1][1], tr3_inv[1][2], tr3_inv[1][3]) ≈ v3ri

    # --- Order 2 (vector6) ---
    ct6 = compile_vector6(fr, 1, 2, 1)
    tr6 = ct6(x)
    v6r = vector6(fr, 1, 2, 1, x)
    @test SVector{3}(tr6[1]...) ≈ SVector{3}(v6r[1], v6r[2], v6r[3])
    @test SVector{3}(tr6[2]...) ≈ SVector{3}(v6r[4], v6r[5], v6r[6])

    # --- Order 3 (vector9) ---
    ct9 = compile_vector9(fr, 1, 2, 1)
    tr9 = ct9(x)
    v9r = vector9(fr, 1, 2, 1, x)
    @test SVector{3}(tr9[1]...) ≈ SVector{3}(v9r[1], v9r[2], v9r[3])
    @test SVector{3}(tr9[2]...) ≈ SVector{3}(v9r[4], v9r[5], v9r[6])
    @test SVector{3}(tr9[3]...) ≈ SVector{3}(v9r[7], v9r[8], v9r[9])

    # --- Order 4 (vector12) ---
    ct12 = compile_vector12(fr, 1, 2, 1)
    tr12 = ct12(x)
    v12r = vector12(fr, 1, 2, 1, x)
    for i in 1:4
        offset = 3*(i-1)
        @test SVector{3}(tr12[i]...) ≈ SVector{3}(v12r[offset+1], v12r[offset+2], v12r[offset+3])
    end

    # Multiple time values
    for t_val in [0.0, 1.0, π/4, 2π]
        tr = compile_vector3(fr, 1, 2, 1)(t_val)
        @test SVector{3}(tr[1][1], tr[1][2], tr[1][3]) ≈ vector3(fr, 1, 2, 1, t_val)
    end
end

@testset "CompiledDirection" begin
    fr = FrameSystem{4, Float64}()
    add_axes!(fr, :ICRF, 1)

    add_point!(fr, :ROOT, 1, 1)

    g(t) = SVector(cos(t), sin(t), 0.0)
    add_point_dynamical!(fr, :Test, 2, 1, 1, g)

    # Add a direction defined in ICRF (axes 1)
    add_direction_position!(fr, :sundir, 1, 2, 1)

    x = π / 3

    # --- Order 1 (direction3) ---
    cd3 = compile_direction3(fr, :sundir)
    tr = cd3(x)
    d3r = direction3(fr, :sundir, 1, x)
    @test SVector{3}(tr[1]...) ≈ d3r

    # --- Order 2 (direction6) ---
    cd6 = compile_direction6(fr, :sundir)
    tr6 = cd6(x)
    d6r = direction6(fr, :sundir, 1, x)
    @test SVector{3}(tr6[1]...) ≈ SVector{3}(d6r[1], d6r[2], d6r[3])
    @test SVector{3}(tr6[2]...) ≈ SVector{3}(d6r[4], d6r[5], d6r[6])

    # --- Error on missing direction ---
    @test_throws ErrorException compile_direction3(fr, :nonexistent)
end

@testset "Compiled non-direct path error" begin
    fr = FrameSystem{2, Float64}()
    add_axes!(fr, :A, 1)
    add_axes_rotating!(fr, :B, 2, 1, t -> angle_to_dcm(t, :Z))
    add_axes_rotating!(fr, :C, 3, 2, t -> angle_to_dcm(t, :X))

    # A → C has path length 3, should error
    @test_throws ErrorException compile_rotation3(fr, 1, 3)
    @test_throws ErrorException compile_rotation6(fr, 1, 3)

    add_point!(fr, :P1, 1, 1)
    f(t) = SVector(cos(t), sin(t), 0.0)
    add_point_dynamical!(fr, :P2, 2, 1, 1, f)
    add_point_dynamical!(fr, :P3, 3, 2, 1, f)

    # P1 → P3 has path length 3, should error
    @test_throws ErrorException compile_vector3(fr, 1, 3, 1)
    @test_throws ErrorException compile_vector6(fr, 1, 3, 1)
end
