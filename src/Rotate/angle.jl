"""
    *(Θ₂::EulerAngles, Θ₁::EulerAngles)

Compute the composed rotation of `Θ₁ -> Θ₂`.

The rotation will be represented by Euler angles (see [`EulerAngles`](@ref))
with the same rotation sequence as `Θ₂`.
"""
@inline function *(Θ₂::EulerAngles, Θ₁::EulerAngles)
    # Convert to quaternions, compute the composition, and convert back to Euler
    # angles.
    q₁ = angle_to_quat(Θ₁)
    q₂ = angle_to_quat(Θ₂)

    return quat_to_angle(q₁ * q₂, Θ₂.seq)
end

"""
    inv(Θ::EulerAngles)

Return the Euler angles that represent the inverse rotation of `Θ`.

The rotation sequence of the result will be the inverse of the input. Hence, if
the input rotation sequence is, for example, `:XYZ`, then the result will be
represented using `:ZYX`.
"""
function inv(Θ::EulerAngles)
    # Check what will be the inverse rotation.
    if Θ.seq == :XYZ
        inv_seq = :ZYX
    elseif Θ.seq == :XZY
        inv_seq = :YZX
    elseif Θ.seq == :YXZ
        inv_seq = :ZXY
    elseif Θ.seq == :YZX
        inv_seq = :XZY
    elseif Θ.seq == :ZXY
        inv_seq = :YXZ
    elseif Θ.seq == :ZYX
        inv_seq = :XYZ
    else
        inv_seq = Θ.seq
    end

    # Return the Euler angle that represented the inverse rotation.
    return EulerAngles(-Θ.a3, -Θ.a2, -Θ.a1, inv_seq)
end

function show(io::IO, Θ::EulerAngles{T}) where T
    # Get if `io` request a compact printing, defaulting to true.
    compact_printing = get(io, :compact, true)

    # Convert the values using `print` and compact printing.
    θ₁_str = sprint(print, Θ.a1; context = :compact => compact_printing)
    θ₂_str = sprint(print, Θ.a2; context = :compact => compact_printing)
    θ₃_str = sprint(print, Θ.a3; context = :compact => compact_printing)
    seq = String(Θ.seq)

    print(io, "EulerAngles{$T}:")
    print(io, " R($seq)  " * θ₁_str * "  " * θ₂_str * "  " * θ₃_str * " rad")

    return nothing
end

Base.eltype(::Type{EulerAngles{T}}) where T = T

