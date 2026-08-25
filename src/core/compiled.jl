
# ==========================================================================================
# Prepared and compiled route APIs. Both resolve graph topology once. Prepared callables
# erase the complete route at one callable boundary; compiled callables retain it in their
# concrete type for maximum local specialization.
# ==========================================================================================

# ------------------------------------------------------------------------------------------
# CompiledRotation
# ------------------------------------------------------------------------------------------

"""
    PreparedRotation{O,T,W}

A compact callable prepared by [`prepare_rotation`](@ref). Its concrete type depends on
the order and scalar signature, but not on route depth or direction.
"""
struct PreparedRotation{O,T,W}
    function_object::W
end

@inline function (rotation::PreparedRotation{O})(t::Number) where {O}
    return Rotation{O}(rotation.function_object(t))
end

"""
    CompiledRotation{O,F}

A fully specialized callable created by [`compile_rotation`](@ref). `F` captures the
complete route, including any inversion, so evaluation can be aggressively inlined.

### Usage
```julia
cr = compile_rotation(fr, :ICRF, :BODY)
R  = cr(t)   # ::Rotation{O} — no FunctionWrapper overhead
```
"""
struct CompiledRotation{O,F}
    function_object::F
end

function CompiledRotation{O}(function_object::F) where {O,F}
    return CompiledRotation{O,F}(function_object)
end

@inline function (rotation::CompiledRotation{O})(t::Number) where {O}
    return Rotation{O}(rotation.function_object(t))
end

# ------------------------------------------------------------------------------------------
# CompiledTranslation
# ------------------------------------------------------------------------------------------

"""
    PreparedTranslation{O,T,W}

A compact callable prepared by [`prepare_translation`](@ref). It returns an
`SVector{3O}` and has a concrete type independent of route depth and direction.
"""
struct PreparedTranslation{O,T,W}
    function_object::W
end

@inline function (translation::PreparedTranslation)(t::Number)
    return translation.function_object(t)
end

"""
    CompiledTranslation{O,F}

A zero-overhead, AD-transparent callable that computes a `Translation{O}` (as `SVector`)
between two points in a specified axes frame. Created via [`compile_translation`](@ref).

### Usage
```julia
ct = compile_translation(fr, :Earth, :Moon, :ICRF)
v  = ct(t)   # ::SVector — no FunctionWrapper overhead
```
"""
struct CompiledTranslation{O,F}
    function_object::F
end

function CompiledTranslation{O}(function_object::F) where {O,F}
    return CompiledTranslation{O,F}(function_object)
end

@inline function (translation::CompiledTranslation{O})(t::Number) where {O}
    return SVector{3 * O}(translation.function_object(t))
end

# ------------------------------------------------------------------------------------------
# CompiledDirection
# ------------------------------------------------------------------------------------------

"""
    PreparedDirection{O,T,W}

A compact callable prepared by [`prepare_direction`](@ref). It returns an
`SVector{3O}` and has a concrete type independent of route depth and direction.
"""
struct PreparedDirection{O,T,W}
    function_object::W
end

@inline function (direction::PreparedDirection)(t::Number)
    return direction.function_object(t)
end

"""
    CompiledDirection{O,F}

A zero-overhead, AD-transparent callable that computes a direction vector at a given time.
Created via [`compile_direction`](@ref).

### Usage
```julia
cd = compile_direction(fr, :SunDir, :ICRF)
d  = cd(t)   # ::SVector — no FunctionWrapper overhead
```
"""
struct CompiledDirection{O,F}
    function_object::F
end

function CompiledDirection{O}(function_object::F) where {O,F}
    return CompiledDirection{O,F}(function_object)
end

@inline function (direction::CompiledDirection{O})(t::Number) where {O}
    return SVector{3 * O}(direction.function_object(t))
end

# A single FunctionWrappersWrapper around each complete operation erases all route details.
function _prepare_rotation(
        rotation::CompiledRotation{O}, ::Type{T}
    ) where {O,T}
    wrapper = _frame_axes_fun_wrapper(Val(O), T, rotation)
    return PreparedRotation{O,T,typeof(wrapper)}(wrapper)
end

function _prepare_translation(
        translation::CompiledTranslation{O}, ::Type{T}
    ) where {O,T}
    wrapper = _frame_vector_fun_wrapper(Val(O), T, translation)
    return PreparedTranslation{O,T,typeof(wrapper)}(wrapper)
end

function _prepare_direction(
        direction::CompiledDirection{O}, ::Type{T}
    ) where {O,T}
    wrapper = _frame_vector_fun_wrapper(Val(O), T, direction)
    return PreparedDirection{O,T,typeof(wrapper)}(wrapper)
end

# ------------------------------------------------------------------------------------------
# compile_rotation
# ------------------------------------------------------------------------------------------

"""
    prepare_rotation(fr::FrameSystem, from, to)
    prepare_rotation(fr::FrameSystem, from, to, ::Val{N})

Resolve an axes route once and return a compact [`PreparedRotation`](@ref). Prepared
rotation types do not encode route depth, direction, or identity. The default order is the
frame system's maximum order; use `Val(N)` to retain only order `N`.

The prepared operation is a snapshot of the graph topology. Prepare it again after changing
the graph.
"""
function prepare_rotation(fr::FrameSystem{O,T}, from, to) where {O,T}
    return prepare_rotation(fr, from, to, Val(O))
end

function prepare_rotation(
        fr::FrameSystem{O,T}, from, to, order::Val{N}
    ) where {O,T,N}
    return _prepare_rotation(compile_rotation(fr, from, to, order), T)
end

"""
    compile_rotation(fr::FrameSystem, from, to)
    compile_rotation(fr::FrameSystem, from, to, ::Val{N})

Compile a zero-overhead rotation callable between axes `from` and `to`. The returned
[`CompiledRotation`](@ref) bypasses `FunctionWrapper` dispatch entirely by extracting
the raw closure and preserving its concrete type.

By default the compiled callable uses the frame system's maximum order `O`. Pass
`Val{N}()` where `N ≤ O` to extract only the `N`-th order closure, avoiding unnecessary
derivative computation in the hot path.

Supports both direct parent-child pairs and multi-hop paths. For multi-hop paths,
the construction loop is type-unstable (closures are composed via `let`), but the
resulting callable is fully typed — so `cr(t)` in a hot loop is zero-overhead.

!!! warning
    The compiled callable captures a snapshot of the current frame graph topology.
    If axes are added after compilation, the callable becomes stale and must be
    recompiled.
"""
function compile_rotation(fr::FrameSystem{O,T}, from, to) where {O,T}
    return compile_rotation(fr, from, to, Val(O))
end

function compile_rotation(
        fr::FrameSystem{O,T}, from, to, ::Val{N}
    ) where {O,T,N}
    N > O && throw(
        ArgumentError("requested order $N exceeds frame system order $O.")
    )

    fromid = axes_id(fr, from)
    toid = axes_id(fr, to)

    if fromid == toid
        return CompiledRotation{N}(t -> Rotation{N}(one(t) * I))
    end

    nodes = _get_axes_nodes(fr, fromid, toid)
    isnothing(nodes) && throw(
        ErrorException("no path between axes $fromid and $toid in the frame system.")
    )

    return _compile_rotation(Val(N), nodes)
end

function _compile_rotation_pair(::Val{N}, from::FrameAxesNode, to::FrameAxesNode) where {N}
    if from.id == to.parentid
        raw_fn = _raw_function(to.f, Val(N))
        return CompiledRotation{N}(raw_fn)
    else
        return let function_object = _raw_function(from.f, Val(N))
            CompiledRotation{N}(
                t -> inv(Rotation{N}(function_object(t))))
        end
    end
end

function _compile_rotation(::Val{N}, nodes::Vector{<:FrameAxesNode}) where {N}
    cr = _compile_rotation_pair(Val(N), nodes[1], nodes[2])

    for i in 3:length(nodes)
        cr_next = _compile_rotation_pair(Val(N), nodes[i-1], nodes[i])
        cr = let inner = cr, outer = cr_next
            CompiledRotation{N}(t -> outer(t) * inner(t))
        end
    end

    return cr
end

# ------------------------------------------------------------------------------------------
# compile_translation
# ------------------------------------------------------------------------------------------

"""
    prepare_translation(fr::FrameSystem, from, to, axes)
    prepare_translation(fr::FrameSystem, from, to, axes, ::Val{N})

Resolve a point and axes route once and return a compact [`PreparedTranslation`](@ref).
The callable returns `SVector{3N}`. Its concrete type does not encode route depth,
direction, or identity. The default order is the frame system's maximum order.

The prepared operation is a snapshot of the graph topology. Prepare it again after changing
the graph.
"""
function prepare_translation(fr::FrameSystem{O,T}, from, to, axes) where {O,T}
    return prepare_translation(fr, from, to, axes, Val(O))
end

function prepare_translation(
        fr::FrameSystem{O,T}, from, to, axes, order::Val{N}
    ) where {O,T,N}
    return _prepare_translation(
        compile_translation(fr, from, to, axes, order), T)
end

"""
    compile_translation(fr::FrameSystem, from, to, axes)
    compile_translation(fr::FrameSystem, from, to, axes, ::Val{N})

Compile a zero-overhead translation callable between points `from` and `to`, expressed in
the given `axes` frame. The returned [`CompiledTranslation`](@ref) bypasses `FunctionWrapper`
dispatch for both the point functions and any axis rotations needed along the path.

By default the compiled callable uses the frame system's maximum order `O`. Pass
`Val{N}()` where `N ≤ O` to extract only the `N`-th order closure.

!!! warning
    The compiled callable captures a snapshot of the current frame graph topology.
    If points or axes are added after compilation, the callable becomes stale and
    must be recompiled.
"""
function compile_translation(fr::FrameSystem{O,T}, from, to, axes) where {O,T}
    return compile_translation(fr, from, to, axes, Val(O))
end

function compile_translation(
        fr::FrameSystem{O,T}, from, to, axes, ::Val{N}
    ) where {O,T,N}
    N > O && throw(
        ArgumentError("requested order $N exceeds frame system order $O.")
    )

    fromid = point_id(fr, from)
    toid = point_id(fr, to)
    axid = axes_id(fr, axes)

    if fromid == toid
        return CompiledTranslation{N}(
            t -> SVector(ntuple(_ -> zero(t), Val(3 * N))))
    end

    nodes = _get_points_nodes(fr, fromid, toid)
    isnothing(nodes) && throw(
        ErrorException("no path between points $fromid and $toid in the frame system.")
    )

    return _compile_translation(Val(N), fr, nodes, axid)
end

function _compile_point_pair(::Val{N}, from::FramePointNode, to::FramePointNode) where {N}
    if from.id == to.parentid
        return to.axesid, _raw_function(to.f, Val(N)), false
    else
        return from.axesid, _raw_function(from.f, Val(N)), true
    end
end

function _compile_translation_hop(::Val{N}, function_object, inverse) where {N}
    inverse && return t -> -Translation{N}(function_object(t))
    return t -> Translation{N}(function_object(t))
end

function _compile_translation(::Val{N}, fr::FrameSystem, nodes::Vector{<:FramePointNode}, axes::Int) where {N}
    if length(nodes) == 2
        axid, raw_fn, inv_flag = _compile_point_pair(Val(N), nodes[1], nodes[2])
        hop = _compile_translation_hop(Val(N), raw_fn, inv_flag)
        if axid == axes
            return CompiledTranslation{N}(let _hop = hop
                t -> SVector(_hop(t))
            end)
        end

        rotation = compile_rotation(fr, axid, axes, Val(N))
        return CompiledTranslation{N}(let _hop = hop, _rotation = rotation
            t -> SVector(_rotation(t) * _hop(t))
        end)
    end

    return _compile_translation_forward(Val(N), fr, nodes, axes)
end

function _compile_translation_forward(::Val{N}, fr::FrameSystem, nodes::Vector{<:FramePointNode}, out_axes::Int) where {N}
    axid1, raw1, inv1 = _compile_point_pair(Val(N), nodes[1], nodes[2])
    first_hop = _compile_translation_hop(Val(N), raw1, inv1)

    compiled_fun = let _hop = first_hop, _axid = axid1
        t -> (_axid, _hop(t))
    end

    prev_axid = axid1

    for i in 3:length(nodes)
        axid_i, raw_i, inv_i = _compile_point_pair(Val(N), nodes[i-1], nodes[i])
        hop_i = _compile_translation_hop(Val(N), raw_i, inv_i)

        if axid_i != prev_axid
            cr = compile_rotation(fr, prev_axid, axid_i, Val(N))
            compiled_fun = let _prev = compiled_fun, _hop = hop_i, _cr = cr, _axid = axid_i
                function (t)
                    _, tr = _prev(t)
                    tr_rotated = _cr(t) * tr
                    return (_axid, tr_rotated + _hop(t))
                end
            end
        else
            compiled_fun = let _prev = compiled_fun, _hop = hop_i, _axid = axid_i
                function (t)
                    _, tr = _prev(t)
                    return (_axid, tr + _hop(t))
                end
            end
        end

        prev_axid = axid_i
    end

    if prev_axid != out_axes
        cr_final = compile_rotation(fr, prev_axid, out_axes, Val(N))
        return CompiledTranslation{N}(let _inner = compiled_fun, _cr = cr_final
            function (t)
                _, tr = _inner(t)
                return SVector(_cr(t) * tr)
            end
        end)
    else
        return CompiledTranslation{N}(let _inner = compiled_fun
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
    prepare_direction(fr::FrameSystem, name::Symbol, axes)
    prepare_direction(fr::FrameSystem, name::Symbol, axes, ::Val{N})

Resolve a direction's axes route once and return a compact [`PreparedDirection`](@ref).
The callable returns `SVector{3N}`. Its concrete type does not encode axes-route depth.
The default order is the frame system's maximum order.

The prepared operation is a snapshot of the graph topology. Prepare it again after changing
the graph.
"""
function prepare_direction(fr::FrameSystem{O,T}, name::Symbol, axes) where {O,T}
    return prepare_direction(fr, name, axes, Val(O))
end

function prepare_direction(
        fr::FrameSystem{O,T}, name::Symbol, axes, order::Val{N}
    ) where {O,T,N}
    return _prepare_direction(compile_direction(fr, name, axes, order), T)
end

"""
    compile_direction(fr::FrameSystem, name::Symbol, axes)
    compile_direction(fr::FrameSystem, name::Symbol, axes, ::Val{N})

Compile a zero-overhead direction callable for the direction `name`, expressed in the
given `axes` frame. The returned [`CompiledDirection`](@ref) bypasses `FunctionWrapper`
dispatch.

By default the compiled callable uses the frame system's maximum order `O`. Pass
`Val{N}()` where `N ≤ O` to extract only the `N`-th order closure.

!!! warning
    The compiled callable captures a snapshot of the current frame graph topology.
    If axes are added after compilation, the callable becomes stale and must be
    recompiled.
"""
function compile_direction(fr::FrameSystem{O,T}, name::Symbol, axes) where {O,T}
    return compile_direction(fr, name, axes, Val(O))
end

function compile_direction(
        fr::FrameSystem{O,T}, name::Symbol, axes, ::Val{N}
    ) where {O,T,N}
    N > O && throw(
        ArgumentError("requested order $N exceeds frame system order $O.")
    )

    if !has_direction(fr, name)
        throw(
            ErrorException("No direction with name $name registered in the frame system.")
        )
    end

    node = directions(fr)[name]
    raw_fn = _raw_function(node.f, Val(N))
    thisaxid = node.axesid
    axid = axes_id(fr, axes)

    if thisaxid != axid
        cr = compile_rotation(fr, thisaxid, axid, Val(N))
        return CompiledDirection{N}(let _fn = raw_fn, _cr = cr
            function (t)
                stv = Translation{N}(_fn(t))
                return SVector(_cr(t) * stv)
            end
        end)
    else
        return CompiledDirection{N}(let _fn = raw_fn
            function (t)
                return SVector(Translation{N}(_fn(t)))
            end
        end)
    end
end
