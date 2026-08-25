# Welcome to FrameTransformations.jl!

_A modern, high-performance and comprehensive set of tools for transformations between any standard and user-defined reference frame._

Are you in search of fundamental routines for efficient and extensible frames transformations?  
If so, this package is the ideal starting point. FrameTransformations.jl is designed to 
provide users with  the ability to create a customized, efficient, flexible, and 
extensible axes/point graph models for mission analysis and space mission design purposes. 

## Features 

- Convert between different time scales and representations (via [Tempo.jl](https://github.com/JuliaSpaceMissionDesign/Tempo.jl)).
- Read binary ephemeris files (via [Ephemerides.jl](https://github.com/JuliaSpaceMissionDesign/Ephemerides.jl) or [CalcephEphemeris.jl](https://github.com/JuliaSpaceMissionDesign/CalcephEphemeris.jl) extensions).
- Create custom reference frame systems with both standard and user-defined points, axes and directions.
- Transform states and their higher-order derivatives between different frames (up to jerk).
- Prepare compact, AD-transparent routes for storage in models via [`prepare_rotation`](@ref), [`prepare_translation`](@ref), and [`prepare_direction`](@ref).
- Compile fully specialized transformations for isolated hot loops via [`compile_rotation`](@ref), [`compile_translation`](@ref), and [`compile_direction`](@ref).

Automatic differentiation is tested through
[DifferentiationInterface.jl](https://github.com/JuliaDiff/DifferentiationInterface.jl)
with ForwardDiff, FiniteDiff, Zygote, and Mooncake backends. Analytic ChainRules rules use
the next available state derivative for direct graph operations. Prepared and compiled
callable rules differentiate their resolved route without traversing the graph.

## Installation 

This package can be installed using Julia's package manager: 
```julia
julia> import Pkg

julia> Pkg.add("FrameTransformations.jl");
```
