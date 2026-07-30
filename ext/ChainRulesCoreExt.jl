module ChainRulesCoreExt

import ChainRulesCore
import ChainRulesCore: rrule
using ChainRulesCore:
    AbstractZero,
    NoTangent,
    ZeroTangent,
    unthunk
import FrameTransformations
using FrameTransformations:
    FrameSystem,
    Rotation,
    direction3,
    direction6,
    direction9,
    direction12,
    order,
    rotation3,
    rotation6,
    rotation9,
    rotation12,
    vector3,
    vector6,
    vector9,
    vector12
using JSMDUtils.Autodiff: derivative
using ReferenceFrameRotations: DCM
using StaticArrays: SMatrix, SVector

# Cotangents arrive as static arrays with Zygote and may arrive as Tangent wrappers
# through other ChainRules consumers such as Mooncake.
to_static_vector(value::SVector) = value
to_static_vector(value::AbstractVector) = SVector(value...)
function to_static_vector(value::ChainRulesCore.Tangent{<:SVector{N}}) where {N}
    return SVector(value.data...)
end
function to_static_vector(
    value::ChainRulesCore.Tangent{<:Any,<:NamedTuple{(:data,)}}
)
    data = value.data
    return data isa ChainRulesCore.Tangent ?
        SVector(data.backing...) : SVector(data...)
end
to_static_vector(value) = SVector(collect(value)...)

dcm_data(value::DCM) = value.data
dcm_data(value::SMatrix{3,3}) = Tuple(value)
dcm_data(value::AbstractMatrix) = Tuple(value)
function dcm_data(value::ChainRulesCore.Tangent)
    data = value.data
    return data isa ChainRulesCore.Tangent ? Tuple(data.backing) : Tuple(data)
end

rotation_components(value::Rotation) = value.m
function rotation_components(value::ChainRulesCore.Tangent)
    components = value.m
    return components isa ChainRulesCore.Tangent ?
        Tuple(components.backing) : Tuple(components)
end
rotation_components(value) = value.m

frobenius_dot(left::NTuple{N}, right::NTuple{N}) where {N} =
    sum(left .* right)
frobenius_dot(left::NTuple, right) = sum(left .* Tuple(right))
frobenius_dot(left, right) = sum(Tuple(left) .* Tuple(right))

contract_output(::AbstractZero, time_derivative) = ZeroTangent()
function contract_output(cotangent, time_derivative::AbstractVector{<:Number})
    values = to_static_vector(unthunk(cotangent))
    return sum(time_derivative .* values)
end
function contract_output(cotangent, time_derivative::Tuple)
    components = rotation_components(unthunk(cotangent))
    result = zero(eltype(dcm_data(first(time_derivative))))
    for index in eachindex(time_derivative)
        result += frobenius_dot(
            dcm_data(time_derivative[index]),
            dcm_data(unthunk(components[index])),
        )
    end
    return result
end
function contract_output(cotangent, time_derivative::AbstractVector{<:DCM})
    return contract_output(cotangent, Tuple(time_derivative))
end

function rule_result(value, time_derivative, ::Val{N}) where {N}
    constant_tangents = ntuple(_ -> NoTangent(), Val(N))

    function frame_pullback(output_tangent)
        time_tangent = contract_output(unthunk(output_tangent), time_derivative)
        return NoTangent(), constant_tangents..., time_tangent
    end

    return value, frame_pullback
end

function fallback_rrule(function_object, arguments...)
    value = function_object(arguments...)
    fixed_arguments = Base.front(arguments)
    time = last(arguments)
    constant_count = length(arguments) - 1

    function fallback_pullback(output_tangent)
        time_derivative = derivative(
            epoch -> function_object(fixed_arguments..., epoch), time
        )
        time_tangent = contract_output(unthunk(output_tangent), time_derivative)
        constant_tangents = ntuple(_ -> NoTangent(), Val(constant_count))
        return NoTangent(), constant_tangents..., time_tangent
    end

    return value, fallback_pullback
end

function vector_rrule(
    function_object,
    higher_order,
    ::Val{N},
    required_order::Int,
    arguments...,
) where {N}
    frame_system = first(arguments)
    if higher_order === nothing || order(frame_system) < required_order
        return fallback_rrule(function_object, arguments...)
    end

    state = higher_order(arguments...)
    value = SVector{N}(ntuple(index -> state[index], Val(N)))
    time_derivative = SVector{N}(
        ntuple(index -> state[index + 3], Val(N))
    )
    return rule_result(value, time_derivative, Val(length(arguments) - 1))
end

function rotation_rrule(
    function_object,
    higher_order,
    ::Val{N},
    required_order::Int,
    arguments...,
) where {N}
    frame_system = first(arguments)
    if higher_order === nothing || order(frame_system) < required_order
        return fallback_rrule(function_object, arguments...)
    end

    state = higher_order(arguments...)
    value = Rotation(ntuple(index -> state[index], Val(N)))
    time_derivative = ntuple(index -> state[index + 1], Val(N))
    return rule_result(value, time_derivative, Val(length(arguments) - 1))
end

for (function_name, higher_order_name, width, required_order) in (
    (:vector3, :vector6, 3, 2),
    (:vector6, :vector9, 6, 3),
    (:vector9, :vector12, 9, 4),
    (:vector12, nothing, 12, 5),
)
    function_object = getfield(FrameTransformations, function_name)
    higher_order = higher_order_name === nothing ?
        nothing : getfield(FrameTransformations, higher_order_name)

    @eval rrule(
        ::typeof($function_object),
        frame_system::FrameSystem,
        from,
        to,
        axes,
        time::Number,
    ) = vector_rrule(
        $function_object,
        $higher_order,
        Val($width),
        $required_order,
        frame_system,
        from,
        to,
        axes,
        time,
    )
end

for (function_name, higher_order_name, width, required_order) in (
    (:direction3, :direction6, 3, 2),
    (:direction6, :direction9, 6, 3),
    (:direction9, :direction12, 9, 4),
    (:direction12, nothing, 12, 5),
)
    function_object = getfield(FrameTransformations, function_name)
    higher_order = higher_order_name === nothing ?
        nothing : getfield(FrameTransformations, higher_order_name)

    @eval rrule(
        ::typeof($function_object),
        frame_system::FrameSystem,
        name::Symbol,
        axes,
        time::Number,
    ) = vector_rrule(
        $function_object,
        $higher_order,
        Val($width),
        $required_order,
        frame_system,
        name,
        axes,
        time,
    )
end

for (function_name, higher_order_name, width, required_order) in (
    (:rotation3, :rotation6, 1, 2),
    (:rotation6, :rotation9, 2, 3),
    (:rotation9, :rotation12, 3, 4),
    (:rotation12, nothing, 4, 5),
)
    function_object = getfield(FrameTransformations, function_name)
    higher_order = higher_order_name === nothing ?
        nothing : getfield(FrameTransformations, higher_order_name)

    @eval rrule(
        ::typeof($function_object),
        frame_system::FrameSystem,
        from,
        to,
        time::Number,
    ) = rotation_rrule(
        $function_object,
        $higher_order,
        Val($width),
        $required_order,
        frame_system,
        from,
        to,
        time,
    )
end

end
