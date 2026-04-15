
# ==========================================================================================
# Compiled fast-path: opt-in API that extracts raw closures from FunctionWrappers,
# giving users zero-overhead, AD-transparent callables for the hot path.
# ==========================================================================================

"""
    _extract_raw_fn(fww::FunctionWrappersWrapper)

Extract the original Julia closure from a `FunctionWrappersWrapper`. The raw closure
preserves its concrete type, enabling the compiler to inline and differentiate through it —
unlike `FunctionWrapper` which erases the type behind a C function pointer.
"""
_extract_raw_fn(fww::FunctionWrappersWrapper) = fww.fw[1].obj[]

# ------------------------------------------------------------------------------------------
# CompiledRotation
# ------------------------------------------------------------------------------------------

"""
    CompiledRotation{O, F}

A zero-overhead, AD-transparent callable that computes a `Rotation{O}` between two axes.
Created via [`compile_rotation`](@ref). The type parameter `F` captures the concrete closure 
type, allowing the compiler to inline the call.

### Usage
```julia
cr = compile_rotation(fr, :ICRF, :BODY)
R  = cr(t)   # ::Rotation{O} — no FunctionWrapper overhead
```
"""
struct CompiledRotation{O,F}
    fun::F
    inverse::Bool
end

function CompiledRotation{O}(fun::F, inverse::Bool) where {O,F}
    return CompiledRotation{O,F}(fun, inverse)
end

function (cr::CompiledRotation{O})(t::Number) where {O}
    R = Rotation{O}(cr.fun(t))
    return cr.inverse ? inv(R) : R
end

# ------------------------------------------------------------------------------------------
# CompiledTranslation
# ------------------------------------------------------------------------------------------

"""
    CompiledTranslation{O, F}

A zero-overhead, AD-transparent callable that computes a `Translation{O}` (as `SVector`)
between two points in a specified axes frame. Created via [`compile_translation`](@ref).

### Usage
```julia
ct = compile_translation(fr, :Earth, :Moon, :ICRF)
v  = ct(t)   # ::SVector — no FunctionWrapper overhead
```
"""
struct CompiledTranslation{O,F}
    fun::F
end

function CompiledTranslation{O}(fun::F) where {O,F}
    return CompiledTranslation{O,F}(fun)
end

function (ct::CompiledTranslation{O})(t::Number) where {O}
    return ct.fun(t)
end

# ------------------------------------------------------------------------------------------
# CompiledDirection
# ------------------------------------------------------------------------------------------

"""
    CompiledDirection{O, F}

A zero-overhead, AD-transparent callable that computes a direction vector at a given time.
Created via [`compile_direction`](@ref).

### Usage
```julia
cd = compile_direction(fr, :SunDir, :ICRF)
d  = cd(t)   # ::SVector — no FunctionWrapper overhead
```
"""
struct CompiledDirection{O,F}
    fun::F
end

function CompiledDirection{O}(fun::F) where {O,F}
    return CompiledDirection{O,F}(fun)
end

function (cd::CompiledDirection{O})(t::Number) where {O}
    return cd.fun(t)
end

# ------------------------------------------------------------------------------------------
# compile_rotation
# ------------------------------------------------------------------------------------------

"""
    compile_rotation(fr::FrameSystem{O}, from, to) where O

Compile a zero-overhead rotation callable between axes `from` and `to`. The returned
[`CompiledRotation`](@ref) bypasses `FunctionWrapper` dispatch entirely by extracting
the raw closure and preserving its concrete type.

Supports both direct parent-child pairs and multi-hop paths. For multi-hop paths,
the construction loop is type-unstable (closures are composed via `let`), but the
resulting callable is fully typed — so `cr(t)` in a hot loop is zero-overhead.

!!! warning
    The compiled callable captures a snapshot of the current frame graph topology.
    If axes are added after compilation, the callable becomes stale and must be
    recompiled.
"""
function compile_rotation(fr::FrameSystem{O,T}, from, to) where {O,T}
    fromid = axes_id(fr, from)
    toid = axes_id(fr, to)

    fromid == toid && return CompiledRotation{O}(_ -> one(T) * I, false)

    nodes = _get_axes_nodes(fr, fromid, toid)
    isnothing(nodes) && throw(
        ErrorException("no path between axes $fromid and $toid in the frame system.")
    )

    return _compile_rotation(Val(O), nodes)
end

function _compile_rotation_pair(::Val{O}, from::FrameAxesNode, to::FrameAxesNode) where {O}
    if from.id == to.parentid
        raw_fn = _extract_raw_fn(to.f[O])
        return CompiledRotation{O}(raw_fn, false)
    else
        raw_fn = _extract_raw_fn(from.f[O])
        return CompiledRotation{O}(raw_fn, true)
    end
end

function _compile_rotation(::Val{O}, nodes::Vector{<:FrameAxesNode}) where {O}
    cr = _compile_rotation_pair(Val(O), nodes[1], nodes[2])

    for i in 3:length(nodes)
        cr_next = _compile_rotation_pair(Val(O), nodes[i-1], nodes[i])
        cr = let inner = cr, outer = cr_next
            CompiledRotation{O}(t -> outer(t) * inner(t), false)
        end
    end

    return cr
end

# ------------------------------------------------------------------------------------------
# compile_translation
# ------------------------------------------------------------------------------------------

"""
    compile_translation(fr::FrameSystem{O}, from, to, axes) where O

Compile a zero-overhead translation callable between points `from` and `to`, expressed in
the given `axes` frame. The returned [`CompiledTranslation`](@ref) bypasses `FunctionWrapper`
dispatch for both the point functions and any axis rotations needed along the path.

!!! warning
    The compiled callable captures a snapshot of the current frame graph topology.
    If points or axes are added after compilation, the callable becomes stale and 
    must be recompiled.
"""
function compile_translation(fr::FrameSystem{O,T}, from, to, axes) where {O,T}
    fromid = point_id(fr, from)
    toid = point_id(fr, to)
    axid = axes_id(fr, axes)

    if fromid == toid
        return CompiledTranslation{O}(_ -> @SVector zeros(T, 3 * O))
    end

    nodes = _get_points_nodes(fr, fromid, toid)
    isnothing(nodes) && throw(
        ErrorException("no path between points $fromid and $toid in the frame system.")
    )

    return _compile_translation(Val(O), fr, nodes, axid)
end

# Extract raw point closure and determine direction (forward/inverse) + axes
function _compile_point_pair(::Val{O}, from::FramePointNode, to::FramePointNode) where {O}
    raw_fn = _extract_raw_fn(to.f[O])
    if from.id == to.parentid
        return to.axesid, raw_fn, false
    else
        raw_fn_from = _extract_raw_fn(from.f[O])
        return from.axesid, raw_fn_from, true
    end
end

function _compile_translation(::Val{O}, fr::FrameSystem, nodes::Vector{<:FramePointNode}, axes::Int) where {O}
    if length(nodes) == 2
        axid, raw_fn, inv_flag = _compile_point_pair(Val(O), nodes[1], nodes[2])
        if axid != axes
            cr = compile_rotation(fr, axid, axes)
            if inv_flag
                return CompiledTranslation{O}(let _fn = raw_fn, _cr = cr
                    t -> SVector(_cr(t) * (-Translation{O}(_fn(t))))
                end)
            else
                return CompiledTranslation{O}(let _fn = raw_fn, _cr = cr
                    t -> SVector(_cr(t) * Translation{O}(_fn(t)))
                end)
            end
        else
            if inv_flag
                return CompiledTranslation{O}(let _fn = raw_fn
                    t -> SVector(-Translation{O}(_fn(t)))
                end)
            else
                return CompiledTranslation{O}(let _fn = raw_fn
                    t -> SVector(Translation{O}(_fn(t)))
                end)
            end
        end
    end

    # Multi-hop: forward pass with compiled rotations where axes change
    return _compile_translation_forward(Val(O), fr, nodes, axes)
end

function _compile_translation_forward(::Val{O}, fr::FrameSystem, nodes::Vector{<:FramePointNode}, out_axes::Int) where {O}
    axid1, raw1, inv1 = _compile_point_pair(Val(O), nodes[1], nodes[2])

    compiled_fun = let _fn = raw1, _inv = inv1
        if _inv
            t -> (axid1, -Translation{O}(_fn(t)))
        else
            t -> (axid1, Translation{O}(_fn(t)))
        end
    end

    prev_axid = axid1

    for i in 3:length(nodes)
        axid_i, raw_i, inv_i = _compile_point_pair(Val(O), nodes[i-1], nodes[i])

        if axid_i != prev_axid
            cr = compile_rotation(fr, prev_axid, axid_i)
            compiled_fun = let _prev = compiled_fun, _fn = raw_i, _inv = inv_i, _cr = cr, _axid = axid_i
                function (t)
                    _, tr = _prev(t)
                    tr_rotated = _cr(t) * tr
                    tr_hop = _inv ? -Translation{O}(_fn(t)) : Translation{O}(_fn(t))
                    return (_axid, tr_rotated + tr_hop)
                end
            end
        else
            compiled_fun = let _prev = compiled_fun, _fn = raw_i, _inv = inv_i, _axid = axid_i
                function (t)
                    _, tr = _prev(t)
                    tr_hop = _inv ? -Translation{O}(_fn(t)) : Translation{O}(_fn(t))
                    return (_axid, tr + tr_hop)
                end
            end
        end

        prev_axid = axid_i
    end

    if prev_axid != out_axes
        cr_final = compile_rotation(fr, prev_axid, out_axes)
        return CompiledTranslation{O}(let _inner = compiled_fun, _cr = cr_final
            function (t)
                _, tr = _inner(t)
                return SVector(_cr(t) * tr)
            end
        end)
    else
        return CompiledTranslation{O}(let _inner = compiled_fun
            function (t)
                _, tr = _inner(t)
                return SVector(tr)
            end
        end)
    end
end

# ------------------------------------------------------------------------------------------
# compile_direction
# ------------------------------------------------------------------------------------------

"""
    compile_direction(fr::FrameSystem{O}, name::Symbol, axes) where O

Compile a zero-overhead direction callable for the direction `name`, expressed in the
given `axes` frame. The returned [`CompiledDirection`](@ref) bypasses `FunctionWrapper`
dispatch.

!!! warning
    The compiled callable captures a snapshot of the current frame graph topology.
    If axes are added after compilation, the callable becomes stale and must be
    recompiled.
"""
function compile_direction(fr::FrameSystem{O}, name::Symbol, axes) where {O}
    if !has_direction(fr, name)
        throw(
            ErrorException("No direction with name $name registered in the frame system.")
        )
    end

    node = directions(fr)[name]
    raw_fn = _extract_raw_fn(node.f[O])
    thisaxid = node.axesid
    axid = axes_id(fr, axes)

    if thisaxid != axid
        cr = compile_rotation(fr, thisaxid, axid)
        return CompiledDirection{O}(let _fn = raw_fn, _cr = cr
            function (t)
                stv = Translation{O}(_fn(t))
                return SVector(_cr(t) * stv)
            end
        end)
    else
        return CompiledDirection{O}(let _fn = raw_fn
            function (t)
                return SVector(Translation{O}(_fn(t)))
            end
        end)
    end
end
