## [Frame System](@id frames_api)

```@docs 
FrameSystem
order 
FrameTransformations.timescale 
```

### Points

```@docs 
has_point
points_graph
points_alias
```

### Axes

```@docs 
has_axes 
axes_graph
axes_alias
```

### Directions

```@docs 
has_direction
directions
```

## [Rotations](@id rotation_api)

```@docs 
Rotation
Base.inv
Translation
```

## [Transformations](@id transformations_api)

### [Points](@id points_transform_api)

```@docs 
vector3
vector6
vector9 
vector12 
```

### [Axes](@id axes_transform_api)

```@docs 
rotation3
rotation6
rotation9
rotation12
```

### [Directions](@id directions_transform_api) 

```@docs 
direction3 
direction6 
direction9
direction12
```

## [Graph, Prepared, and Compiled Operations](@id compiled_api)

| Use case | API |
| --- | --- |
| One-off graph query | `rotation6`, `vector6`, `direction6` |
| Store in a model or trajectory | `prepare_rotation`, `prepare_translation`, `prepare_direction` |
| Maximum local throughput | `compile_rotation`, `compile_translation`, `compile_direction` |

Direct operations retain graph flexibility. Prepared operations resolve topology once and
place the complete route behind one compact callable boundary, keeping large SciML model
types independent of route depth. Compiled operations retain the complete route in their
concrete type so Julia can maximize inlining, at the cost of more downstream specialization
and compilation.

Preparation and compilation snapshot graph topology. Additions to the graph do not change an
existing callable; prepare or compile it again after changing the graph. The default form uses
the frame system's maximum order. Pass `Val(N)` to store only order `N`.

Rotation callables return `Rotation{N}`. Translation and direction callables return
`SVector{3N}`. Scalar types promote with generic numerical times, including ForwardDiff duals.
Differentiate direct, prepared, and compiled callables through DifferentiationInterface with
an `AutoForwardDiff()` backend.

### Types

```@docs 
PreparedRotation
PreparedTranslation
PreparedDirection
CompiledRotation
CompiledTranslation
CompiledDirection
```

### Constructors

```@docs 
prepare_rotation
prepare_translation
prepare_direction
compile_rotation
compile_translation
compile_direction
```
