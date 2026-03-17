"""
    CompiledRotation{O, F}

A pre-compiled, allocation-free rotation function extracted from a FrameSystem.
Bypasses FunctionWrappers and graph path lookups by storing the raw closure directly.

Type parameter `O` is the rotation order and `F` is the closure type.
"""
struct CompiledRotation{O, F}
    fun::F
    inverse::Bool
end

@inline function (cr::CompiledRotation{O})(t::Number) where O
    R = cr.fun(t)
    return cr.inverse ? inv(R) : R
end

"""
    CompiledTranslation{O, F}

A pre-compiled, allocation-free translation function extracted from a FrameSystem.
Bypasses FunctionWrappers and graph path lookups by storing the raw closure directly.

Type parameter `O` is the translation order and `F` is the closure type.
"""
struct CompiledTranslation{O, F}
    fun::F
    axesid::Int
    inverse::Bool
end

@inline function (ct::CompiledTranslation{O})(t::Number) where O
    tr = ct.fun(t)
    return ct.inverse ? -tr : tr
end

"""
    CompiledDirection{O, F}

A pre-compiled, allocation-free direction function extracted from a FrameSystem.
Bypasses FunctionWrappers and Dict lookups by storing the raw closure directly.

Type parameter `O` is the direction order and `F` is the closure type.
"""
struct CompiledDirection{O, F}
    fun::F
    axesid::Int
end

@inline function (cd::CompiledDirection{O})(t::Number) where O
    return Translation{O}(cd.fun(t))
end

# ==========================================================================================
# Rotations
# ==========================================================================================

"""
    compile_rotation(fr::FrameSystem, from, to, ::Val{order}) -> CompiledRotation

Extract the raw rotation closure for a direct parent-child axes pair.
Only works for directly connected axes (path length = 2).

The `order` parameter selects which derivative level to compile (1-4):
- `Val(1)` → position only (rotation3)
- `Val(2)` → position + velocity (rotation6)
- `Val(3)` → position + velocity + acceleration (rotation9)
- `Val(4)` → position + velocity + acceleration + jerk (rotation12)
"""
function compile_rotation(fr::FrameSystem{O,T}, from, to, ::Val{order}) where {O,T,order}
    fromid = axes_id(fr, from)
    toid = axes_id(fr, to)

    path = get_path(axes_graph(fr), fromid, toid)
    length(path) == 2 || error(
        "compile_rotation only supports direct parent-child pairs (path length must be 2, got $(length(path)))"
    )

    f1 = get_mappednode(axes_graph(fr), path[1])
    f2 = get_mappednode(axes_graph(fr), path[2])

    if f1.id == f2.parentid
        raw_fn = f2.f.fun[order].fw[1].obj[]
        return CompiledRotation{order, typeof(raw_fn)}(raw_fn, false)
    else
        raw_fn = f1.f.fun[order].fw[1].obj[]
        return CompiledRotation{order, typeof(raw_fn)}(raw_fn, true)
    end
end

compile_rotation3(fr::FrameSystem, from, to)  = compile_rotation(fr, from, to, Val(1))
compile_rotation6(fr::FrameSystem, from, to)  = compile_rotation(fr, from, to, Val(2))
compile_rotation9(fr::FrameSystem, from, to)  = compile_rotation(fr, from, to, Val(3))
compile_rotation12(fr::FrameSystem, from, to) = compile_rotation(fr, from, to, Val(4))

# ==========================================================================================
# Translations (vectors / point positions)
# ==========================================================================================

"""
    compile_translation(fr::FrameSystem, from_point, to_point, axes, ::Val{order}) -> CompiledTranslation

Extract the raw translation closure for a direct parent-child point pair.
Only works for directly connected points (path length = 2).

The `order` parameter selects which derivative level to compile (1-4):
- `Val(1)` → position only (vector3)
- `Val(2)` → position + velocity (vector6)
- `Val(3)` → position + velocity + acceleration (vector9)
- `Val(4)` → position + velocity + acceleration + jerk (vector12)
"""
function compile_translation(fr::FrameSystem{O,T}, from_point, to_point, axes, ::Val{order}) where {O,T,order}
    fromid = point_id(fr, from_point)
    toid = point_id(fr, to_point)

    path = get_path(points_graph(fr), fromid, toid)
    length(path) == 2 || error(
        "compile_translation only supports direct parent-child pairs (path length must be 2, got $(length(path)))"
    )

    p1 = get_mappednode(points_graph(fr), path[1])
    p2 = get_mappednode(points_graph(fr), path[2])

    axid = axes_id(fr, axes)

    if p1.id == p2.parentid
        raw_fn = p2.f.fun[order].fw[1].obj[]
        return CompiledTranslation{order, typeof(raw_fn)}(raw_fn, p2.axesid, false)
    else
        raw_fn = p1.f.fun[order].fw[1].obj[]
        return CompiledTranslation{order, typeof(raw_fn)}(raw_fn, p1.axesid, true)
    end
end

compile_vector3(fr::FrameSystem, from, to, axes)  = compile_translation(fr, from, to, axes, Val(1))
compile_vector6(fr::FrameSystem, from, to, axes)  = compile_translation(fr, from, to, axes, Val(2))
compile_vector9(fr::FrameSystem, from, to, axes)  = compile_translation(fr, from, to, axes, Val(3))
compile_vector12(fr::FrameSystem, from, to, axes) = compile_translation(fr, from, to, axes, Val(4))

# ==========================================================================================
# Directions
# ==========================================================================================

"""
    compile_direction(fr::FrameSystem, name::Symbol, ::Val{order}) -> CompiledDirection

Extract the raw direction closure from a FrameSystem.

The `order` parameter selects which derivative level to compile (1-4):
- `Val(1)` → direction only (direction3)
- `Val(2)` → direction + rate (direction6)
- `Val(3)` → direction + rate + acceleration (direction9)
- `Val(4)` → direction + rate + acceleration + jerk (direction12)

!!! note
    Unlike rotations and translations, directions are stored in a Dict rather than a graph,
    so there is no path-length constraint.  However, if the output axes differ from the
    definition axes, a rotation is still needed at runtime.  The compiled direction only
    bypasses the Dict lookup and FunctionWrapper call — the caller is responsible for
    any axes rotation.
"""
function compile_direction(fr::FrameSystem{O,T}, name::Symbol, ::Val{order}) where {O,T,order}
    if !has_direction(fr, name)
        error("No direction with name :$name registered in the frame system.")
    end
    node = directions(fr)[name]
    raw_fn = node.f.fun[order].fw[1].obj[]
    return CompiledDirection{order, typeof(raw_fn)}(raw_fn, node.axesid)
end

compile_direction3(fr::FrameSystem, name::Symbol)  = compile_direction(fr, name, Val(1))
compile_direction6(fr::FrameSystem, name::Symbol)  = compile_direction(fr, name, Val(2))
compile_direction9(fr::FrameSystem, name::Symbol)  = compile_direction(fr, name, Val(3))
compile_direction12(fr::FrameSystem, name::Symbol) = compile_direction(fr, name, Val(4))
