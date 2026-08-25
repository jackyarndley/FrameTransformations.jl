using FrameTransformations
using DifferentiationInterface: AutoForwardDiff, derivative, second_derivative
using ReferenceFrameRotations
using StaticArrays
import ForwardDiff
using Test

const backend = AutoForwardDiff()

@testset "FrameSystem" begin 
    fr = FrameSystem{4, Float64}()
    add_axes!(fr, :ICRF, 1)

    r(t) = angle_to_dcm(t, :Z)
    dr(t) = derivative(r, backend, t)
    d2r(t) = second_derivative(r, backend, t)
    d3r(t) = derivative(d2r, backend, t)
    add_axes_rotating!(fr, :Test, 2, 1, r)

    f(t) = SVector(cos(t), sin(t), 1.0)
    df(t) = derivative(f, backend, t)
    d2f(t) = second_derivative(f, backend, t)
    d3f(t) = derivative(d2f, backend, t)

    add_point!(fr, :ROOT, 1, 1)
    add_point_dynamical!(fr, :Test, 2, 1, 1, f)
    
    x = π/3
    @test vector3(fr, 1, 2, 1, x) == f(x)
    @test vector6(fr, 1, 2, 1, x) == vcat( f(x), df(x) )
    @test vector9(fr, 1, 2, 1, x) == vcat( f(x), df(x), d2f(x) )
    @test vector12(fr, 1, 2, 1, x) == vcat( f(x), df(x), d2f(x), d3f(x) )

    x = π/3
    @test rotation3(fr, 1, 2, x)[1] == r(x)
    @test rotation6(fr, 1, 2, x)[2] == dr(x)
    @test rotation9(fr, 1, 2, x)[3] == d2r(x)
    @test rotation12(fr, 1, 2, x)[4] == d3r(x)
end
