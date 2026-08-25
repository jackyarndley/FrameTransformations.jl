using DifferentiationInterface: AutoForwardDiff, derivative
using FrameTransformations
using IERSConventions
using ReferenceFrameRotations
using Tempo
using Test
import ForwardDiff

const backend = AutoForwardDiff()
const time = 123_456_789.0

@testset "Orekit TEME definition" begin
    centuries = time / Tempo.CENTURY2SEC
    correction = IERSConventions.eop_δΔψ(iers2010b, centuries)
    equation = IERSConventions.equation_equinoxes(
        iers2010b, centuries, correction)
    expected = angle_to_dcm(equation, :Z) *
               iers_rot3_gcrf_to_tod(time, iers2010b)

    @test teme_rot3_gcrf_to_teme(time) ≈ expected atol=2e-15 rtol=2e-15

    j2000_reference = DCM([
        0.9999999996369622 -6.056843713980243e-8 2.6945724773411362e-5
        5.981382710408946e-8 0.9999999996078635 2.800481486491666e-5
        -2.6945726459052826e-5 -2.8004813243018796e-5 0.9999999992448292
    ])
    @test teme_rot3_gcrf_to_teme(0.0) ≈ j2000_reference atol=1e-14 rtol=1e-14

    simplified_equation = first(iers_nutation_comp(iers2010b, centuries)) *
                          cos(iers_obliquity(iers2010b, centuries))
    @test abs(equation - simplified_equation) > 1e-12
end

@testset "TEME cumulative orders and registration" begin
    rotation_functions = (
        teme_rot3_gcrf_to_teme,
        teme_rot6_gcrf_to_teme,
        teme_rot9_gcrf_to_teme,
        teme_rot12_gcrf_to_teme,
    )
    for order in 1:4
        value = rotation_functions[order](time)
        rotation = order == 1 ? (value,) : value
        @test length(rotation) == order
        @test all(component -> component isa DCM, rotation)
    end

    rotation6_value = teme_rot6_gcrf_to_teme(time)
    @test derivative(teme_rot3_gcrf_to_teme, backend, time) ≈
          rotation6_value[2] atol=2e-14 rtol=2e-12

    frames = FrameSystem{4,Float64,TerrestrialTime}()
    add_axes_gcrf!(frames)
    @test_nowarn add_axes_teme!(frames, :TEME, AXESID_GCRF)
    for order in 1:4
        direct = (rotation3, rotation6, rotation9, rotation12)[order](
            frames, :GCRF, :TEME, time)
        value = rotation_functions[order](time)
        expected = order == 1 ? (value,) : value
        @test all(isapprox.(direct.m, expected; atol=2e-14, rtol=2e-12))
        prepared = prepare_rotation(frames, :GCRF, :TEME, Val(order))(time)
        @test all(isapprox.(prepared.m, expected; atol=2e-14, rtol=2e-12))
    end

    @test_throws ArgumentError add_axes_teme!(
        FrameSystem{1,Float64,TerrestrialTime}(), :TEME, 999)
end

@testset "TEME prepares the Tempo conversion once" begin
    frames = FrameSystem{2,Float64}()
    add_axes_gcrf!(frames)
    add_axes_teme!(frames, :TEME, AXESID_GCRF)

    to_tt = Tempo.prepare_time_conversion(TDB, TT)
    @test rotation3(frames, :GCRF, :TEME, time)[1] ≈
          teme_rot3_gcrf_to_teme(to_tt(time)) atol=2e-15 rtol=2e-15

    prepared = prepare_rotation(frames, :GCRF, :TEME, Val(2))
    @test derivative(epoch -> prepared(epoch)[1], backend, time) ≈
          prepared(time)[2] atol=2e-14 rtol=2e-12
end
