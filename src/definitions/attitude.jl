
"""
    add_axes_fixed_quaternion!(frames::FrameSystem, name::Symbol, id::Int, parent, q::Quaternion)
   
Add axes `name` with id `id` to `frames` with a fixed-offset from `parent`. 
Fixed offset axes have a constant orientation with respect to their `parent` axes, 
represented by the quaternion `q`.

### See also 
See also [`add_axes_fixedoffset!`](@ref).
"""
function add_axes_fixed_quaternion!(
    frames::FrameSystem, name::Symbol, id::Int, parent, q::Quaternion
)
    add_axes_fixedoffset!(frames, name, id, parent, quat_to_dcm(q))
end

"""
    add_axes_fixed_angles!(frames, name::Symbol, id::Int, parent,
        angles::AbstractVector{N}, seq::Symbol)
   
Add axes `name` with id `id` to `frames` with a fixed-offset from `parent`. 
Fixed offset axes have a constant orientation with respect to their `parent` axes, 
represented by Euler `angles`.

The rotation sequence is defined by `seq` specifing the rotation axes. The possible
values depends on the number of rotations as follows:

- **1 rotation** (`θ₁`): `:X`, `:Y`, or `:Z`.
- **2 rotations** (`θ₁`, `θ₂`): `:XY`, `:XZ`, `:YX`, `:YZ`, `:ZX`, or `:ZY`.
- **3 rotations** (`θ₁`, `θ₂`, `θ₃`): `:XYX`, `XYZ`, `:XZX`, `:XZY`, `:YXY`, `:YXZ`, `:YZX`,
    `:YZY`, `:ZXY`, `:ZXZ`, `:ZYX`, or `:ZYZ`

### See also 
See also [`add_axes_fixedoffset!`](@ref).
"""
function add_axes_fixed_angles!(
    frames::FrameSystem,
    name::Symbol,
    id::Int,
    parent,
    angles::AbstractVector{N},
    seq::Symbol,
) where {N}
    add_axes_fixedoffset!(
        frames, name, id, parent, angle_to_dcm(angles..., seq)
    )
end

"""
    add_axes_fixed_angleaxis!(frames, name::Symbol, id::Int, parent,
        angle::Number, axis::AbstractVector{N})
   
Add axes `name` with id `id` to `frames` with a fixed-offset from `parent`. 
Fixed offset axes have a constant orientation with respect to their `parent` axes, 
represented by Euler `angle` [rad] and Euler `axis`.

### See also 
See also [`add_axes_fixedoffset!`](@ref).
"""
function add_axes_fixed_angleaxis!(
    frames::FrameSystem,
    name::Symbol,
    id::Int,
    parent,
    angle::Number,
    axis::AbstractVector{N},
) where {N}
    normalized_axis = unitvec(axis)
    add_axes_fixedoffset!(
        frames, name, id, parent, angleaxis_to_dcm(angle, normalized_axis)
    )
end
