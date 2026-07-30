
"""
    add_direction!(frames, name::Symbol, axes, fun,
        first_derivative=nothing, second_derivative=nothing, third_derivative=nothing)

Add a new direction node to `frames`. The orientation of these direction depends only 
on time and is computed through the custom functions provided by the user. 

The input functions must accept only time as argument and their outputs must be as follows: 

- `fun`: return a direction vector.
- `first_derivative`: return a direction vector and its first time derivative.
- `second_derivative`: return a direction vector and its first two time derivatives.
- `third_derivative`: return a direction vector and its first three time derivatives.

Missing derivative functions are computed via automatic differentiation.

!!! warning 
    It is expected that the input functions and their outputs have the correct signature. This 
    function does not perform any checks on the output types. 
"""
function add_direction!(
    frames::FrameSystem{O,N}, name::Symbol, axes, fun,
    first_derivative=nothing, second_derivative=nothing, third_derivative=nothing
) where {O,N}
    has_direction(frames, name) && throw(
        ArgumentError("direction with name=$name is already registered in the frame system.")
    )

    for (derivative_order, function_object) in enumerate(
        (first_derivative, second_derivative, third_derivative)
    )
        if O < derivative_order + 1 && !isnothing(function_object)
            @warn "ignoring $function_object, frame system order is less than $(derivative_order + 1)"
        end
    end

    funs = DirectionFunctions{O,N}(
        t -> Translation{O}(fun(t)),

        # First derivative
        if isnothing(first_derivative)
            t -> _translation_from_parts(
                Val(O), fun(t), derivative1(fun, t)
            )
        else
            t -> Translation{O}(first_derivative(t))
        end,

        # Second derivative
        if isnothing(second_derivative)
            (
                if isnothing(first_derivative)
                    t -> _translation_from_parts(
                        Val(O),
                        fun(t),
                        derivative1(fun, t),
                        derivative2(fun, t),
                    )
                else
                    t -> Translation{O}(
                        SVector(
                            first_derivative(t)...,
                            derivative2(fun, t)...,
                        )
                    )
                end
            )
        else
            t -> Translation{O}(second_derivative(t))
        end,

        # Third derivative 
        if isnothing(third_derivative)
            (
                if isnothing(second_derivative)
                    (
                        if isnothing(first_derivative)
                            t -> _translation_from_parts(
                                Val(O),
                                fun(t),
                                derivative1(fun, t),
                                derivative2(fun, t),
                                derivative3(fun, t),
                            )
                        else
                            t -> Translation{O}(
                                SVector(
                                    first_derivative(t)...,
                                    derivative2(fun, t)...,
                                    derivative3(fun, t)...,
                                )
                            )
                        end
                    )
                else
                    t -> Translation{O}(
                        SVector(
                            second_derivative(t)...,
                            derivative3(fun, t)...,
                        )
                    )
                end
            )
        else
            t -> Translation{O}(third_derivative(t))
        end,
    )

    axid = axes_id(frames, axes)
    has_axes(frames, axid) || throw(
        ArgumentError("Axes with ID $axid are not registered in the input frame system")
    )
    dir = DirectionDefinition{O,N}(name, length(directions(frames)) + 1, axid, funs)
    push!(directions(frames), Pair(name, dir))
    nothing
end
