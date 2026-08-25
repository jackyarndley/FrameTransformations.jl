
# FrameTransformations.jl

_A modern high-performance set of tools for transformations between standard and user-defined reference frame._

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaspacemissiondesign.github.io/FrameTransformations.jl/stable/) 
[![Dev Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaspacemissiondesign.github.io/FrameTransformations.jl/dev/) 
[![Build Status](https://github.com/JuliaSpaceMissionDesign/FrameTransformations.jl/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/JuliaSpaceMissionDesign/FrameTransformations.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/JuliaSpaceMissionDesign/FrameTransformations.jl/branch/main/graph/badge.svg?token=7fj9BjJhKF)](https://codecov.io/gh/JuliaSpaceMissionDesign/FrameTransformations.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

Are you in search of fundamental routines for efficient and extensible frames transformations?  
If so, this package is the ideal starting point. FrameTransformations.jl is designed to 
provide users with  the ability to create a customized, efficient, flexible, and 
extensible axes/point graph models for mission analysis and space mission design purposes. 

## Features 

- Convert between different time scales and representations (via [Tempo.jl](https://github.com/JuliaSpaceMissionDesign/Tempo.jl));
- Read binary ephemeris files (via [Ephemerides.jl](https://github.com/JuliaSpaceMissionDesign/Ephemerides.jl) or [CalcephEphemeris.jl](https://github.com/JuliaSpaceMissionDesign/CalcephEphemeris.jl))
- Create custom reference frame systems with both standard and user-defined points, axes and directions.
- Transform states and their higher-order derivatives between different frames (up to jerk).
- Prepare compact, AD-transparent routes for storage in models via `prepare_rotation`, `prepare_translation`, and `prepare_direction`.
- Compile fully specialized transformations for isolated hot loops via `compile_rotation`, `compile_translation`, and `compile_direction`.

Automatic differentiation is tested through
[DifferentiationInterface.jl](https://github.com/JuliaDiff/DifferentiationInterface.jl)
with its ForwardDiff backend. Direct, prepared, and compiled callables preserve dual-number
times and differentiate through their resolved numerical operations.

## Transformation modes

| Use case | API |
| --- | --- |
| One-off graph query | `rotation6`, `vector6`, `direction6` |
| Store in a model or trajectory | `prepare_rotation`, `prepare_translation`, `prepare_direction` |
| Maximum local throughput | `compile_rotation`, `compile_translation`, `compile_direction` |

Prepared and compiled operations snapshot graph topology, so graph changes require preparing
or compiling again. Use `Val(N)` to select an exact kinematic order. Rotations return
`Rotation{N}`; translations and directions return `SVector{3N}`. Prepared types erase route
depth to limit downstream SciML specialization, while compiled types expose the complete route
to Julia's optimizer.

## Installation 

This package can be installed using Julia's package manager: 
```julia
julia> import Pkg

julia> Pkg.add("FrameTransformations.jl");
```

## Documentation 
For further information on this package and its tutorials please refer to the 
[stable documentation](https://juliaspacemissiondesign.github.io/FrameTransformations.jl/stable/).

## Benchmarking

The repository includes simple benchmark and profiling scripts under `benchmark/`.
Run the benchmark suite from the repository root with:

```julia
julia benchmark/runbenchmarks.jl
```

Measure package load, frame construction, route preparation/compilation, first evaluation,
and first composite-model use in a clean Julia process with:

```julia
julia --startup-file=no benchmark/runfirstuse.jl
```

To run the SPICE position/state workloads used by
[Brahe PR #376](https://github.com/duncaneddy/brahe/pull/376), use:

```julia
julia benchmark/runspicebenchmarks.jl
```

This compares prepared `Ephemerides.jl` chains with direct, prepared, and compiled
`FrameTransformations.jl` paths for single Sun/Moon/Mars-barycenter queries, a
10,000-epoch sequential propagation pattern, and a Sun + Moon third-body acceleration
adapter. Kernel loading, route preparation, frame-graph construction, and transformation
compilation are excluded from the measurements. The bundled DE432s test kernel is used by
default; set `SPICE_KERNEL` to a DE440s BSP file for the exact kernel used by Brahe:

```bash
SPICE_KERNEL=/path/to/de440s.bsp julia benchmark/runspicebenchmarks.jl
```

The benchmark script uses `BenchmarkTools.jl` and compares direct, prepared, and compiled paths
for synthetic multi-hop examples and the DE440 lunar-frame rotation case at frame-system
orders 2, 3, and 4. The DE440 suite includes both rotation and Earth-to-Moon vector
benchmarks expressed in the lunar ME421 frame, plus DifferentiationInterface evaluations.
The benchmark and test environments are tracked as Julia workspace projects from the root
`Project.toml` on Julia 1.12+, with the legacy path-based setup retained for older Julia
releases.

You can shorten or lengthen the run with:

```bash
BENCHMARK_SECONDS=1 BENCHMARK_SAMPLES=10 julia benchmark/runbenchmarks.jl
```

For a quick CPU profile of representative direct, prepared, and compiled paths, run:

```julia
julia benchmark/profile.jl
```

## Support
If you found this package useful, please consider starring the repository. We also encourage 
you to take a look at other astrodynamical packages of the [JSMD](https://github.com/JuliaSpaceMissionDesign/) organisation.
