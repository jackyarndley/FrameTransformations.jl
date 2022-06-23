export δdcm, orthonormalize

# StaticArray.jl API
# ==============================================================================

@inline Base.@propagate_inbounds function getindex(dcm::DCM, i::Int)
    return dcm.data[i]
end

function Tuple(dcm::DCM)
    return dcm.data
end

function similar_type(::Type{A}, ::Type{T}, s::Size{(3,3)}) where {A<:DCM, T}
    return DCM{T}
end

# Julia API
# ==============================================================================

function summary(io::IO, ::DCM{T}) where T
    print(io, "DCM{" * string(T) * "}")
end

"""
    orthonormalize(dcm::DCM)

Perform the Gram-Schmidt orthonormalization process in the `dcm` and return the
new matrix.

!!! warning
    This function does not check if the columns of the input matrix span a
    three-dimensional space. If not, then the returned matrix should have `NaN`.
    Notice, however, that such input matrix is not a valid direction cosine
    matrix.
"""
function orthonormalize(dcm::DCM)
    e₁ = dcm[:, 1]
    e₂ = dcm[:, 2]
    e₃ = dcm[:, 3]

    en₁  = e₁ / norm(e₁)
    enj₂ = e₂ - (en₁ ⋅ e₂) * en₁
    en₂  = enj₂ / norm(enj₂)
    enj₃ = e₃ - (en₁ ⋅ e₃) * en₁
    enj₃ = enj₃ - (en₂ ⋅ enj₃) * en₂
    en₃  = enj₃ / norm(enj₃)

    return DCM(hcat(en₁, en₂, en₃))
end

#### Kinematics

"""
    δdcm(Dba::DCM, wba_b::AbstractArray)

Compute the time-derivative of the `dcm` that rotates a reference frame `a` into
alignment with the reference frame `b` in which the angular velocity of `b` with
respect to `a`, and represented in `b`, is `wba_b`.
"""
function δdcm(Dba::DCM, wba_b::AbstractArray)
    # Auxiliary variable.
    w = wba_b

    # Check the dimensions.
    if length(wba_b) != 3
        throw(ArgumentError("The angular velocity vector must have three components."))
    end

    wx = SMatrix{3, 3}(
          0  , -w[3], +w[2],
        +w[3],   0  , -w[1],
        -w[2], +w[1],   0,
    )'

    # Return the time-derivative.
    return -wx * Dba
end
