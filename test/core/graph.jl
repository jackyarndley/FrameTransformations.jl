using FrameTransformations

using Tempo
using LinearAlgebra
using ReferenceFrameRotations: DCM
using StaticArrays
using Test 

@testset "Point node" begin
    begin
        fps = FrameTransformations.FramePointFunctions{1, Float64}()
        fp = FrameTransformations.FramePointNode{1, Float64}(:P0, 1, 0, 1, fps)

        @test fp.name == :P0
        @test fp.id == 1 
        @test fp.parentid == 0
        @test fp.axesid == 1
        @test fp.f === fps
    end

    begin 
        f1(t) = Translation{2}(SA[cos(t), sin(t), 1.0])
        f2(t) = Translation{2}(SA[cos(t), sin(t), 1.0, -sin(t), cos(t), 0.0])
        fps = FrameTransformations.FramePointFunctions{2, Float64}(f1, f2)

        @test Translation{1}(f1(π/3)) == fps[1].fw[1](π/3)
        @test f2(π/3) == fps[2].fw[1](π/3)
        @test f1(0)[2] == zeros(3)
        @test fps[Val(1)](0.1) isa Translation{1}
        @test fps[Val(2)](0.1) isa Translation{2}
        @test FrameTransformations._raw_function(fps, Val(1))(0.1) isa Translation{1}
        @test FrameTransformations._raw_function(fps, Val(2))(0.1) isa Translation{2}
    end
end;

@testset "Exact-order function storage" begin
    point_calls = zeros(Int, 4)
    point_functions = ntuple(Val(4)) do index
        t -> begin
            point_calls[index] += 1
            Translation{index,Float64}()
        end
    end
    points = FrameTransformations.FramePointFunctions{4,Float64}(point_functions...)

    axes_calls = zeros(Int, 4)
    axes_functions = ntuple(Val(4)) do index
        t -> begin
            axes_calls[index] += 1
            Rotation{index}(one(t) * I)
        end
    end
    axes = FrameTransformations.FrameAxesFunctions{4,Float64}(axes_functions...)

    for index in 1:4
        @test order(points[Val(index)](0.0)) == index
        @test order(axes[Val(index)](0.0)) == index
        @test order(FrameTransformations._raw_function(points, Val(index))(0.0)) == index
        @test order(FrameTransformations._raw_function(axes, Val(index))(0.0)) == index
    end
    @test point_calls == fill(2, 4)
    @test axes_calls == fill(2, 4)

    # A lower-order query touches only its corresponding provider.
    fill!(point_calls, 0)
    fill!(axes_calls, 0)
    points[Val(1)](0.0)
    axes[Val(1)](0.0)
    @test point_calls == [1, 0, 0, 0]
    @test axes_calls == [1, 0, 0, 0]
end;

@testset "Axes node" begin
    begin 
        faxs = FrameTransformations.FrameAxesFunctions{1, Float64}()
        fax = FrameTransformations.FrameAxesNode{1, Float64}(:Ax0, 1, 1, faxs)

        @test fax.name == :Ax0 
        @test fax.id == 1
        @test fax.parentid == 1
        @test fax.f === faxs
        @test faxs[1](0.0f0) == faxs[1](0.0)
    end
end;

@testset "FrameSystem" begin 
    g = FrameSystem{2, Float64}()

    @test order(g) == 2 
    @test FrameTransformations.timescale(g) == BarycentricDynamicalTime
    @test points_graph(g) === g.points.graph
    @test axes_graph(g) === g.axes.graph
    @test points_alias(g) === g.points.alias
    @test axes_alias(g) === g.axes.alias

    @test !has_point(g, 1)
    @test !has_axes(g, 1)
end;

@testset "Lazy route caches" begin
    frames = FrameSystem{1,Float64}()
    add_axes!(frames, :A, 1)
    add_axes_fixedoffset!(frames, :B, 2, :A, DCM{Float64}(I))
    add_point!(frames, :P, 1, :A)
    add_point_fixedoffset!(frames, :Q, 2, :P, :A, SA[1.0, 0.0, 0.0])

    @test isempty(frames.axes_route_cache)
    @test isempty(frames.points_route_cache)

    rotation3(frames, :A, :B, 0.0)
    @test collect(keys(frames.axes_route_cache)) == [(1, 2)]

    vector3(frames, :P, :Q, :A, 0.0)
    @test collect(keys(frames.points_route_cache)) == [(1, 2)]

    add_axes_fixedoffset!(frames, :C, 3, :B, DCM{Float64}(I))
    @test isempty(frames.axes_route_cache)
    add_point_fixedoffset!(frames, :R, 3, :Q, :A, SA[0.0, 1.0, 0.0])
    @test isempty(frames.points_route_cache)
end;
