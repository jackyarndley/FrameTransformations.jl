
export skew 

"""
    skew(a)

Create a skew matrix from the first three elements of `a`.
"""
skew(a) = SMatrix{3, 3}(0, a[3], -a[2], -a[3], 0, a[1], a[2], -a[1], 0)


"""
    angle_to_δdcm(θ::Number[, ϕ::Number[, γ::Number]], rot_seq::Symbol = :ZYX)

Compute the derivative of the direction cosine matrix that perform a set of rotations (`θ`, 
`ϕ`, `γ`) about the coordinate axes in `rot_seq`. 

The rotation sequence is defined by a `Symbol` specifing the rotation axes. The possible 
values depend on the number of rotations as follows: 

- **1 rotation** (`θ₁`): `:X`, `:Y`, or `:Z`.
- **2 rotations** (`θ₁`, `θ₂`): `:XY`, `:XZ`, `:YX`, `:YZ`, `:ZX`, or `:ZY`.
- **3 rotations** (`θ₁`, `θ₂`, `θ₃`): `:XYX`, `XYZ`, `:XZX`, `:XZY`, `:YXY`,
    `:YXZ`, `:YZX`, `:YZY`, `:ZXY`, `:ZXZ`, `:ZYX`, or `:ZYZ`

!!! note 
    This function assigns `dcm = A3 * A2 * A1` in which `Ai` is the DCM related with
    the *i*-th rotation, `i Є [1,2,3]`. If the *i*-th rotation is not specified,
    then `Ai = I`.

### See also 
See also [`angle_to_δ²dcm`](@ref)
"""
function angle_to_δdcm(θ, rot_seq::Symbol)
    
    s, c = sincos(θ[1]).*θ[2]

    if rot_seq == :X 
        return DCM(0, 0, 0, 0, -s, -c, 0, c, -s)

    elseif rot_seq == :Y 
        return DCM(-s, 0, c, 0, 0, 0, -c, 0, -s)

    elseif rot_seq == :Z
        return DCM(-s, -c, 0, c, -s, 0, 0, 0, 0)

    else
        throw(ArgumentError("rot_seq must be :X, :Y, or :Z"))
    end

end

function angle_to_δdcm(θ, ϕ, rot_seq::Symbol)

    s, c = sincos(θ[1])
    b, a = sincos(ϕ[1])

    δϕ, δθ = ϕ[2], θ[2]

    if rot_seq == :XY
        return DCM(
            -b*δϕ, 0, a*δϕ, 
            a*s*δϕ + b*c*δθ, -s*δθ, b*s*δϕ-a*c*δθ, 
            b*s*δθ -a*c*δϕ, c*δθ, -b*c*δϕ-a*s*δθ
        )
        
    elseif rot_seq == :XZ
        return DCM(
            -b*δϕ, -a*δϕ, 0, 
            a*c*δϕ -b*s*δθ, -b*c*δϕ-a*s*δθ, -c*δθ, 
            a*s*δϕ+b*c*δθ, a*c*δθ-b*s*δϕ, -s*δθ
        )

    elseif rot_seq == :YX 
        return DCM(
            -δθ*s, δθ*c*b + δϕ*s*a, δθ*c*a - δϕ*s*b, 
            0, -δϕ*b, -δϕ*a, 
            -δθ*c, δϕ*c*a - δθ*s*b, -δϕ*b*c - δθ*a*s
        )

    elseif rot_seq == :YZ 
        return DCM(
            -a*s*δθ-b*c*δϕ, b*s*δθ-a*c*δϕ, c*δθ, 
            a*δϕ, -b*δϕ, 0, 
            b*s*δϕ-a*c*δθ, b*c*δθ+a*s*δϕ, -s*δθ
        )

    elseif rot_seq == :ZX 
        return DCM(
            -s*δθ, b*s*δϕ-a*c*δθ, b*c*δθ+a*s*δϕ, 
            c*δθ, -b*c*δϕ-a*s*δθ, b*s*δθ-a*c*δϕ, 
            0, a*δϕ, -b*δϕ        
        )

    elseif rot_seq == :ZY 
        return DCM(
            -a*s*δθ-b*c*δϕ, -c*δθ, a*c*δϕ-b*s*δθ, 
            a*c*δθ -b*s*δϕ, -s*δθ, b*c*δθ+a*s*δϕ, 
            -a*δϕ, 0, -b*δϕ
        )
    
    else 
        throw(ArgumentError(
            "The rotation sequence :$rot_seq is not valid."
        ))
    end

end

function angle_to_δdcm(θ, ϕ, γ, rot_seq::Symbol)

    s, c = sincos(θ[1])
    b, a = sincos(ϕ[1])
    e, d = sincos(γ[1])

    δγ, δϕ, δθ = γ[2], ϕ[2], θ[2]

    if rot_seq == :ZYX
        return DCM(
            -δθ*a*s         - δϕ*b*c, 
            -δθ*(d*c+e*b*s) + δϕ*e*a*c  + δγ*(e*s+d*b*c),
             δθ*(e*c-d*b*s) + δϕ*d*a*c  + δγ*(d*s-e*b*c),
             δθ*a*c         - δϕ*b*s,
             δθ*(e*b*c-d*s) + δϕ*e*a*s  + δγ*(d*b*s-e*c),
             δθ*(e*s+d*b*c) + δϕ*a*d*s  - δγ*(d*c+e*b*s),
                            - δϕ*a, 
                            - δϕ*e*b    + δγ*d*a, 
                            - δϕ*d*b    - δγ*e*a
         )

    elseif rot_seq == :XYX
        return DCM(
                              -δϕ*b, 
                              δϕ*e*a    + δγ*d*b, 
                              δϕ*d*a    - δγ*e*b, 
             δθ*b*c         + δϕ*a*s, 
            -δθ*(d*s+e*a*c) + δϕ*e*b*s  - δγ*(e*c+d*a*s),
             δθ*(e*s-d*a*c) + δϕ*d*b*s  + δγ*(e*a*s-d*c),
             δθ*b*s         - δϕ*a*c, 
             δθ*(d*c-e*a*s) - δϕ*e*b*c  + δγ*(d*a*c-e*s),
            -δθ*(e*c+d*a*s) - δϕ*d*b*c  - δγ*(d*s+e*a*c)
        )

    elseif rot_seq == :XYZ
        return DCM(
                            - δϕ*d*b    - δγ*e*a, 
                              δϕ*e*b    - δγ*d*a, 
                              δϕ*a, 
             δθ*(d*b*c-e*s) + δϕ*d*a*s  + δγ*(d*c-e*b*s),
            -δθ*(e*b*c+d*s) - δϕ*e*a*s  - δγ*(d*b*s+e*c),
            -δθ*a*c         + δϕ*b*s,
             δθ*(d*b*s+e*c) - δϕ*d*a*c  + δγ*(e*b*c+d*s),
             δθ*(d*c-e*b*s) + δϕ*e*a*c  + δγ*(d*b*c-e*s),
            -δθ*a*s         - δϕ*b*c
        )

    elseif rot_seq == :XZX
        return DCM(
                            - δϕ*b, 
                            - δϕ*d*a    + δγ*e*b, 
                              δϕ*e*a    + δγ*d*b, 
            -δθ*b*s         + δϕ*a*c,
            -δθ*(d*a*s+e*c) - δϕ*d*b*c  - δγ*(e*a*c+d*s),
             δθ*(e*a*s-d*c) + δϕ*e*b*c  + δγ*(e*s-d*a*c),
             δθ*b*c         + δϕ*a*s, 
             δθ*(d*a*c-e*s) - δϕ*d*b*s  + δγ*(d*c-e*a*s),
            -δθ*(e*a*c+d*s) + δϕ*e*b*s  - δγ*(d*a*s+e*c)
        )

    elseif rot_seq == :XZY
        return DCM(
                            - δϕ*d*b    - δγ*e*a, 
                            - δϕ*a,  
                            - δϕ*e*b    + δγ*d*a, 
             δθ*(e*c-d*b*s) + δϕ*d*a*c  + δγ*(d*s-e*b*c), 
            -δθ*a*s         - δϕ*b*c,  
            -δθ*(d*c+e*b*s) + δϕ*e*a*c  + δγ*(d*b*c+e*s),
             δθ*(d*b*c+e*s) + δϕ*d*a*s  - δγ*(e*b*s+d*c),
             δθ*a*c         - δϕ*b*s, 
             δθ*(e*b*c-d*s) + δϕ*e*a*s  + δγ*(d*b*s-e*c)
        )

    elseif rot_seq == :YXY
        return DCM(
            -δθ*(d*s+e*a*c) + δϕ*e*b*s  - δγ*(e*c+d*a*s),
             δθ*b*c         + δϕ*a*s,
             δθ*(d*a*c-e*s) - δϕ*d*b*s  + δγ*(d*c-e*a*s),
                              δϕ*e*a    + δγ*d*b,
                            - δϕ*b,
                            - δϕ*d*a    + δγ*e*b,
             δθ*(e*a*s-d*c) + δϕ*e*b*c  + δγ*(e*s-d*a*c),
            -δθ*b*s         + δϕ*a*c,
            -δθ*(e*c+d*a*s) - δϕ*d*b*c  - δγ*(d*s+e*a*c)

        )

    elseif rot_seq == :YXZ
        return DCM(
             δθ*(e*b*c-d*s) + δϕ*e*a*s  + δγ*(d*b*s-e*c),
             δθ*(e*s+d*b*c) + δϕ*d*a*s  - δγ*(d*c+e*b*s),
             δθ*a*c         - δϕ*b*s, 
                            - δϕ*e*b    + δγ*d*a,
                            - δϕ*d*b    - δγ*e*a,
                            - δϕ*a, 
            -δθ*(d*c+e*b*s) + δϕ*e*a*c  + δγ*(e*s+d*b*c),
             δθ*(e*c-d*b*s) + δϕ*d*a*c  + δγ*(d*s-e*b*c),
            -δθ*a*s         - δϕ*b*c, 
        )

    elseif rot_seq == :YZX
        return DCM(
            -δθ*a*s         - δϕ*b*c, 
             δθ*(d*b*s+e*c) - δϕ*d*a*c  + δγ*(e*b*c+d*s),
             δθ*(d*c-e*b*s) + δϕ*e*a*c  + δγ*(d*b*c-e*s),
                              δϕ*a,
                            - δϕ*d*b    - δγ*e*a,
                              δϕ*e*b    - δγ*d*a,
            -δθ*a*c         + δϕ*b*s,
             δθ*(d*b*c-e*s) + δϕ*d*a*s  + δγ*(d*c-e*b*s),
            -δθ*(e*b*c+d*s) - δϕ*e*a*s  - δγ*(d*b*s+e*c)
        )

    elseif rot_seq == :YZY
        return DCM(
            -δθ*(e*c+d*a*s) - δϕ*d*b*c  - δγ*(e*a*c+d*s),
             δθ*b*s         - δϕ*a*c, 
             δθ*(d*c-e*a*s) - δϕ*e*b*c  + δγ*(d*a*c-e*s),
                              δϕ*d*a    - δγ*e*b,
                            - δϕ*b, 
                              δϕ*e*a    + δγ*d*b,
             δθ*(e*s-d*a*c) + δϕ*d*b*s  + δγ*(e*a*s-d*c),
             δθ*b*c         + δϕ*a*s, 
            -δθ*(e*a*c+d*s) + δϕ*e*b*s  - δγ*(d*a*s+e*c)
        )

    elseif rot_seq == :ZXY
        return DCM(
            -δθ*(d*s+e*b*c) - δϕ*e*a*s  - δγ*(e*c+d*b*s),
            -δθ*a*c         + δϕ*b*s,
             δθ*(d*b*c-e*s) + δϕ*d*a*s  + δγ*(d*c-e*b*s),
             δθ*(d*c-e*b*s) + δϕ*e*a*c  + δγ*(d*b*c-e*s),
            -δθ*a*s         - δϕ*b*c,
             δθ*(e*c+d*b*s) - δϕ*d*a*c  + δγ*(d*s+e*b*c),
                              δϕ*e*b    - δγ*d*a,
                              δϕ*a,
                            - δϕ*d*b    - δγ*e*a
        )

    elseif rot_seq == :ZXZ
        return DCM(
            -δθ*(d*s+e*a*c) + δϕ*e*b*s  - δγ*(e*c+d*a*s),
             δθ*(e*s-d*a*c) + δϕ*d*b*s  + δγ*(e*a*s-d*c),
             δθ*b*c         + δϕ*a*s,
             δθ*(d*c-e*a*s) - δϕ*e*b*c  + δγ*(d*a*c-e*s),
            -δθ*(e*c+d*a*s) - δϕ*d*b*c  - δγ*(d*s+e*a*c),
             δθ*b*s         - δϕ*a*c,
                              δϕ*a*e    + δγ*d*b,
                              δϕ*d*a    - δγ*e*b,
                            - δϕ*b
        )
        
    elseif rot_seq == :ZYZ
        return DCM(
            -δθ*(d*a*s+e*c) - δϕ*d*b*c  - δγ*(e*a*c+d*s),
             δθ*(e*a*s-d*c) + δϕ*e*b*c  + δγ*(e*s-d*a*c),
            -δθ*b*s         + δϕ*a*c,
             δθ*(d*a*c-e*s) - δϕ*d*b*s  + δγ*(d*c-e*a*s),
            -δθ*(e*a*c+d*s) + δϕ*e*b*s  - δγ*(d*a*s+e*c),
             δθ*b*c         + δϕ*a*s,
                            - δϕ*a*d    + δγ*e*b,
                              δϕ*e*a    + δγ*d*b,
                            - δϕ*b
        )

    else
        throw(ArgumentError("The rotation sequence :$rot_seq is not valid."))
    end
end


"""
    angle_to_δ²dcm(θ::Number[, ϕ::Number[, γ::Number]], rot_seq::Symbol = :ZYX)

Compute the second order time derivative of the direction cosine matrix that perform 
a set of rotations (`θ`, `ϕ`, `γ`) about the coordinate axes in `rot_seq`. 

The rotation sequence is defined by a `Symbol` specifing the rotation axes. The possible 
values depend on the number of rotations as follows: 

- **1 rotation** (`θ₁`): `:X`, `:Y`, or `:Z`.
- **2 rotations** (`θ₁`, `θ₂`): `:XY`, `:XZ`, `:YX`, `:YZ`, `:ZX`, or `:ZY`.
- **3 rotations** (`θ₁`, `θ₂`, `θ₃`): `:XYX`, `XYZ`, `:XZX`, `:XZY`, `:YXY`,
    `:YXZ`, `:YZX`, `:YZY`, `:ZXY`, `:ZXZ`, `:ZYX`, or `:ZYZ`

!!! note 
    This function assigns `dcm = A3 * A2 * A1` in which `Ai` is the DCM related with
    the *i*-th rotation, `i Є [1,2,3]`. If the *i*-th rotation is not specified,
    then `Ai = I`.

### See also 
See also [`angle_to_δdcm`](@ref)
"""
function angle_to_δ²dcm(θ, rot_seq::Symbol)

    s, c = sincos(θ[1])

    a, b = (s, c).*θ[2]^2
    d, e = (s, c).*θ[3]

    if rot_seq == :X 
        return DCM(0, 0, 0, 0, -d-b, a-e, 0, e-a, -d-b )

    elseif rot_seq == :Y 
        return DCM(-d-b, 0, e-a, 0, 0, 0, a-e, 0, -d-b)

    elseif rot_seq == :Z
        return DCM(-d-b, a-e, 0, e-a, -d-b,  0, 0, 0, 0) 

    else
        throw(ArgumentError("rot_seq must be :X, :Y, or :Z"))
    end

end


function angle_to_δ²dcm(θ, ϕ, rot_seq::Symbol)

    s, c = sincos(θ[1])
    b, a = sincos(ϕ[1])

    δϕ, δ²ϕ = ϕ[2], ϕ[3]
    δθ, δ²θ = θ[2], θ[3]

    δθ², δϕ² = δθ^2, δϕ^2

    δϕθ = 2*δϕ*δθ
    A = (δθ²+ δϕ²)

    if rot_seq == :XY
        return DCM(
                        - δ²ϕ*b     - δϕ²*a, 
            0,
                          δ²ϕ*a     - δϕ²*b,
             δ²θ*c*b    + δ²ϕ*a*s   - A*b*s     + δϕθ*a*c,
            -δ²θ*s                  - δθ²*c,
            -δ²θ*a*c    + δ²ϕ*b*s   + A*a*s     + δϕθ*c*b,
             δ²θ*b*s    - δ²ϕ*a*c   + A*c*b     + δϕθ*a*s,
             δ²θ*c                  - δθ²*s,
            -δ²θ*a*s    - δ²ϕ*c*b   - A*a*c     + δϕθ*b*s
        )

    elseif rot_seq == :XZ
        return DCM(
                        - δ²ϕ*b     - δϕ²*a,
                        - δ²ϕ*a     + δϕ²*b,
            0,
            -δ²θ*b*s    + δ²ϕ*a*c   - A*c*b     - δϕθ*a*s,
            -δ²θ*a*s    - δ²ϕ*c*b   - A*a*c     + δϕθ*b*s,
            -δ²θ*c                  + δθ²*s ,
             δ²θ*c*b    + δ²ϕ*a*s   - A*b*s     + δϕθ*a*c,
             δ²θ*a*c    - δ²ϕ*b*s   - A*a*s     - δϕθ*c*b,
            -δ²θ*s                  - δθ²*c
        )
        
        
    elseif rot_seq == :YX 
        return DCM(
            -δ²θ*s                  - δθ²*c,
             δ²θ*c*b    + δ²ϕ*a*s   - A*b*s     + δϕθ*a*c,
             δ²θ*a*c    - δ²ϕ*b*s   - A*a*s     - δϕθ*c*b,
            0,
                        - δ²ϕ*b     - δϕ²*a,
                        - δ²ϕ*a     + δϕ²*b,
            -δ²θ*c                  + δθ²*s,
            -δ²θ*b*s    + δ²ϕ*a*c   - A*c*b     - δϕθ*a*s,
            -δ²θ*a*s    - δ²ϕ*c*b   - A*a*c     + δϕθ*b*s 
        )

    elseif rot_seq == :YZ 
        return DCM(
            -δ²θ*a*s    - δ²ϕ*c*b   - A*a*c     + δϕθ*b*s,
             δ²θ*b*s    - δ²ϕ*a*c   + A*c*b     + δϕθ*a*s,
             δ²θ*c                  - δθ²*s,
                          δ²ϕ*a     - δϕ²*b,
                        - δ²ϕ*b     - δϕ²*a ,
            0,
            -δ²θ*a*c    + δ²ϕ*b*s   + A*a*s     + δϕθ*c*b,
            +δ²θ*c*b    + δ²ϕ*a*s   - A*b*s     + δϕθ*a*c,
            -δ²θ*s                  - δθ²*c
        )

    elseif rot_seq == :ZX 
        return DCM(
            -δ²θ*s                  - δθ²*c,
            -δ²θ*a*c    + δ²ϕ*b*s   + A*a*s     + δϕθ*c*b,
             δ²θ*c*b    + δ²ϕ*a*s   - A*b*s     + δϕθ*a*c,
             δ²θ*c                  - δθ²*s,
            -δ²θ*a*s    - δ²ϕ*c*b   - A*a*c     + δϕθ*b*s,
             δ²θ*b*s    - δ²ϕ*a*c   + A*c*b     + δϕθ*a*s,
            0,
                          δ²ϕ*a     - δϕ²*b,
                        - δ²ϕ*b     - δϕ²*a
        )

    elseif rot_seq == :ZY 
        return DCM(
            -δ²θ*a*s    - δ²ϕ*c*b   - A*a*c     + δϕθ*b*s,
            -δ²θ*c                  + δθ²*s,
            -δ²θ*b*s    + δ²ϕ*a*c   - A*c*b     - δϕθ*a*s,
             δ²θ*a*c    - δ²ϕ*b*s   - A*a*s     - δϕθ*c*b,
            -δ²θ*s                  - δθ²*c,
             δ²θ*c*b    + δ²ϕ*a*s   - A*b*s + δϕθ*a*c,
                        - δ²ϕ*a     + δϕ²*b,
            0,
                        - δ²ϕ*b     - δϕ²*a 
        )
    
    else 
        throw(ArgumentError(
            "The rotation sequence :$rot_seq is not valid."
        ))
    end

end


function angle_to_δ²dcm(θ, ϕ, γ, rot_seq::Symbol)

    s, c = sincos(θ[1])
    b, a = sincos(ϕ[1])
    e, d = sincos(γ[1])

    δθ, δ²θ = θ[2], θ[3]
    δϕ, δ²ϕ = ϕ[2], ϕ[3]
    δγ, δ²γ = γ[2], γ[3]

    δθ², δϕ², δγ² = δθ^2, δϕ^2, δγ^2

    δϕθ = 2*δϕ*δθ
    δϕγ = 2*δϕ*δγ
    δθγ = 2*δθ*δγ

    A = (δθ² + δϕ² + δγ²)
    B = δθ² + δϕ²
    C = δθ² + δγ²
    D = δγ² + δϕ²

	if rot_seq == :ZYX
		return DCM(
			-a*(c*B + s*δ²θ) - b*(c*δ²ϕ - s*δϕθ),
            c*(d*(a*δϕγ + δ²γ*b - δ²θ) + e*(a*δ²ϕ - A*b + δθγ)) + s*(e*(δ²γ - a*δϕθ - δ²θ*b) + d*(C-δθγ*b)),
            c*(a*(d*δ²ϕ - e*δϕγ) + e*(δ²θ - δ²γ*b) + d*(δθγ- A*b)) + s*(d*(δ²γ-a*δϕθ -δ²θ*b) - e*(C-δθγ*b)),
			a*(c*δ²θ - B*s) - b*(c*δϕθ + s*δ²ϕ),
            c*(e*(a*δϕθ -δ²γ + b*δ²θ) + d*(δθγ*b-C)) + s*(d*(a*δϕγ + δ²γ*b - δ²θ) +e*(a*δ²ϕ + δθγ - A*b)),
            c*d*(a*δϕθ-δ²γ+δ²θ*b) + d*s*(a*δ²ϕ- A*b + δθγ) -e*s*(a*δϕγ - δ²θ + δ²γ*b) + c*e*(C-δθγ*b),
			-a*δ²ϕ + b*δϕ²,
            a*(d*δ²γ-D*e) - b*(d*δϕγ + e*δ²ϕ),
			-a*(D*d + e*δ²γ) + b*(-d*δ²ϕ + e*δϕγ),
		)
	elseif rot_seq == :XYX
		return DCM(
			-a*δϕ² - b*δ²ϕ,
			a*(d*δϕγ + e*δ²ϕ) + b*(d*δ²γ - D*e),
			a*(d*δ²ϕ - e*δϕγ) - b*(D*d + e*δ²γ),
			a*(c*δϕθ + s*δ²ϕ) + b*(c*δ²θ - B*s),
            e*(s*(A*a + b*δ²ϕ + δθγ) + c*(b*δϕθ - δ²γ - δ²θ*a)) + d*(s*(b*δϕγ - δ²γ*a - δ²θ) -c*(C + δθγ*a)),
            d*(s*(A*a+b*δ²ϕ+δθγ) + c*(b*δϕθ-δ²γ-a*δ²θ)) + e*(s*(δ²γ*a- b*δϕγ+δ²θ) + c*(C+δθγ*a)),
			a*(s*δϕθ-c*δ²ϕ ) + b*(B*c + s*δ²θ),
            c*(e*(-A*a-b*δ²ϕ-δθγ) + d*(δ²θ+δ²γ*a- b*δϕγ)) + s*(e*(b*δϕθ-δ²γ-a*δ²θ) - d*(C+a*δθγ)),
            c*(d*(-A*a-δθγ) + b*(-d*δ²ϕ+e*δϕγ ) - e*(δ²θ+δ²γ*a)) + s*(d*(b*δϕθ-δ²γ- δ²θ*a) + e*(C+δθγ*a)),
		)
	elseif rot_seq == :XYZ # TODO: da qui in giu raggruppa in fattori
		return DCM(
            e*(b*δϕγ-a*δ²γ) - d*(D*a+b*δ²ϕ),
            +a*(e*D-d*δ²γ) + b*(d*δϕγ + e*δ²ϕ),
			a*δ²ϕ - b*δϕ²,
			a*c*d*δϕθ + a*d*s*δ²ϕ - a*e*s*δϕγ - b*d*s*δϕ² + δ²γ*(-b*e*s + c*d) + δ²θ*(b*c*d - e*s) + δγ²*(-b*d*s - c*e) + δθ²*(-b*d*s - c*e) + δθγ*(-b*c*e - d*s),
			-a*c*e*δϕθ - a*d*s*δϕγ - a*e*s*δ²ϕ + b*e*s*δϕ² + δ²γ*(-b*d*s - c*e) + δ²θ*(-b*c*e - d*s) + δγ²*(b*e*s - c*d) + δθ²*(b*e*s - c*d) + δθγ*(-b*c*d + e*s),
			-a*(c*δ²θ - s*B) + b*(c*δϕθ + s*δ²ϕ),
			-a*c*d*δ²ϕ + a*c*e*δϕγ + a*d*s*δϕθ + b*c*d*δϕ² + δ²γ*(b*c*e + d*s) + δ²θ*(b*d*s + c*e) + δγ²*(b*c*d - e*s) + δθ²*(b*c*d - e*s) + δθγ*(-b*e*s + c*d),
			a*c*d*δϕγ + a*c*e*δ²ϕ - a*e*s*δϕθ - b*c*e*δϕ² + δ²γ*(b*c*d - e*s) + δ²θ*(-b*e*s + c*d) + δγ²*(-b*c*e - d*s) + δθ²*(-b*c*e - d*s) + δθγ*(-b*d*s - c*e),
			-a*(c*B + s*δ²θ) - b*(c*δ²ϕ - s*δϕθ),
		)
	elseif rot_seq == :XZX
		return DCM(
			-a*δϕ² - b*δ²ϕ,
			-a*d*δ²ϕ + a*e*δϕγ + b*d*δγ² + b*d*δϕ² + b*e*δ²γ,
			a*d*δϕγ + a*e*δ²ϕ + b*d*δ²γ - b*e*δγ² - b*e*δϕ²,
			a*c*δ²ϕ - a*s*δϕθ - b*c*δθ² - b*c*δϕ² - b*s*δ²θ,
			-a*c*d*δϕ² - b*c*d*δ²ϕ + b*c*e*δϕγ + b*d*s*δϕθ + δ²γ*(-a*c*e - d*s) + δ²θ*(-a*d*s - c*e) + δγ²*(-a*c*d + e*s) + δθ²*(-a*c*d + e*s) + δθγ*(a*e*s - c*d),
			a*c*e*δϕ² + b*c*d*δϕγ + b*c*e*δ²ϕ - b*e*s*δϕθ + δ²γ*(-a*c*d + e*s) + δ²θ*(a*e*s - c*d) + δγ²*(a*c*e + d*s) + δθ²*(a*c*e + d*s) + δθγ*(a*d*s + c*e),
			a*c*δϕθ + a*s*δ²ϕ + b*c*δ²θ - b*s*δθ² - b*s*δϕ²,
			-a*d*s*δϕ² - b*c*d*δϕθ - b*d*s*δ²ϕ + b*e*s*δϕγ + δ²γ*(-a*e*s + c*d) + δ²θ*(a*c*d - e*s) + δγ²*(-a*d*s - c*e) + δθ²*(-a*d*s - c*e) + δθγ*(-a*c*e - d*s),
			a*e*s*δϕ² + b*c*e*δϕθ + b*d*s*δϕγ + b*e*s*δ²ϕ + δ²γ*(-a*d*s - c*e) + δ²θ*(-a*c*e - d*s) + δγ²*(a*e*s - c*d) + δθ²*(a*e*s - c*d) + δθγ*(-a*c*d + e*s),
		)
	elseif rot_seq == :XZY
		return DCM(
			-a*d*δγ² - a*d*δϕ² - a*e*δ²γ - b*d*δ²ϕ + b*e*δϕγ,
			-a*δ²ϕ + b*δϕ²,
			a*d*δ²γ - a*e*δγ² - a*e*δϕ² - b*d*δϕγ - b*e*δ²ϕ,
			a*c*d*δ²ϕ - a*c*e*δϕγ - a*d*s*δϕθ - b*c*d*δϕ² + δ²γ*(-b*c*e + d*s) + δ²θ*(-b*d*s + c*e) + δγ²*(-b*c*d - e*s) + δθ²*(-b*c*d - e*s) + δθγ*(b*e*s + c*d),
			-a*c*δθ² - a*c*δϕ² - a*s*δ²θ - b*c*δ²ϕ + b*s*δϕθ,
			a*c*d*δϕγ + a*c*e*δ²ϕ - a*e*s*δϕθ - b*c*e*δϕ² + δ²γ*(b*c*d + e*s) + δ²θ*(-b*e*s - c*d) + δγ²*(-b*c*e + d*s) + δθ²*(-b*c*e + d*s) + δθγ*(-b*d*s + c*e),
			a*c*d*δϕθ + a*d*s*δ²ϕ - a*e*s*δϕγ - b*d*s*δϕ² + δ²γ*(-b*e*s - c*d) + δ²θ*(b*c*d + e*s) + δγ²*(-b*d*s + c*e) + δθ²*(-b*d*s + c*e) + δθγ*(-b*c*e + d*s),
			a*c*δ²θ - a*s*δθ² - a*s*δϕ² - b*c*δϕθ - b*s*δ²ϕ,
			a*c*e*δϕθ + a*d*s*δϕγ + a*e*s*δ²ϕ - b*e*s*δϕ² + δ²γ*(b*d*s - c*e) + δ²θ*(b*c*e - d*s) + δγ²*(-b*e*s - c*d) + δθ²*(-b*e*s - c*d) + δθγ*(b*c*d + e*s),
		)
	elseif rot_seq == :YXY
		return DCM(
			a*e*s*δϕ² + b*c*e*δϕθ + b*d*s*δϕγ + b*e*s*δ²ϕ + δ²γ*(-a*d*s - c*e) + δ²θ*(-a*c*e - d*s) + δγ²*(a*e*s - c*d) + δθ²*(a*e*s - c*d) + δθγ*(-a*c*d + e*s),
			a*c*δϕθ + a*s*δ²ϕ + b*c*δ²θ - b*s*δθ² - b*s*δϕ²,
			-a*d*s*δϕ² - b*c*d*δϕθ - b*d*s*δ²ϕ + b*e*s*δϕγ + δ²γ*(-a*e*s + c*d) + δ²θ*(a*c*d - e*s) + δγ²*(-a*d*s - c*e) + δθ²*(-a*d*s - c*e) + δθγ*(-a*c*e - d*s),
			a*d*δϕγ + a*e*δ²ϕ + b*d*δ²γ - b*e*δγ² - b*e*δϕ²,
			-a*δϕ² - b*δ²ϕ,
			-a*d*δ²ϕ + a*e*δϕγ + b*d*δγ² + b*d*δϕ² + b*e*δ²γ,
			a*c*e*δϕ² + b*c*d*δϕγ + b*c*e*δ²ϕ - b*e*s*δϕθ + δ²γ*(-a*c*d + e*s) + δ²θ*(a*e*s - c*d) + δγ²*(a*c*e + d*s) + δθ²*(a*c*e + d*s) + δθγ*(a*d*s + c*e),
			a*c*δ²ϕ - a*s*δϕθ - b*c*δθ² - b*c*δϕ² - b*s*δ²θ,
			-a*c*d*δϕ² - b*c*d*δ²ϕ + b*c*e*δϕγ + b*d*s*δϕθ + δ²γ*(-a*c*e - d*s) + δ²θ*(-a*d*s - c*e) + δγ²*(-a*c*d + e*s) + δθ²*(-a*c*d + e*s) + δθγ*(a*e*s - c*d),
		)
	elseif rot_seq == :YXZ
		return DCM(
			a*c*e*δϕθ + a*d*s*δϕγ + a*e*s*δ²ϕ - b*e*s*δϕ² + δ²γ*(b*d*s - c*e) + δ²θ*(b*c*e - d*s) + δγ²*(-b*e*s - c*d) + δθ²*(-b*e*s - c*d) + δθγ*(b*c*d + e*s),
			a*c*d*δϕθ + a*d*s*δ²ϕ - a*e*s*δϕγ - b*d*s*δϕ² + δ²γ*(-b*e*s - c*d) + δ²θ*(b*c*d + e*s) + δγ²*(-b*d*s + c*e) + δθ²*(-b*d*s + c*e) + δθγ*(-b*c*e + d*s),
			a*c*δ²θ - a*s*δθ² - a*s*δϕ² - b*c*δϕθ - b*s*δ²ϕ,
			a*d*δ²γ - a*e*δγ² - a*e*δϕ² - b*d*δϕγ - b*e*δ²ϕ,
			-a*d*δγ² - a*d*δϕ² - a*e*δ²γ - b*d*δ²ϕ + b*e*δϕγ,
			-a*δ²ϕ + b*δϕ²,
			a*c*d*δϕγ + a*c*e*δ²ϕ - a*e*s*δϕθ - b*c*e*δϕ² + δ²γ*(b*c*d + e*s) + δ²θ*(-b*e*s - c*d) + δγ²*(-b*c*e + d*s) + δθ²*(-b*c*e + d*s) + δθγ*(-b*d*s + c*e),
			a*c*d*δ²ϕ - a*c*e*δϕγ - a*d*s*δϕθ - b*c*d*δϕ² + δ²γ*(-b*c*e + d*s) + δ²θ*(-b*d*s + c*e) + δγ²*(-b*c*d - e*s) + δθ²*(-b*c*d - e*s) + δθγ*(b*e*s + c*d),
			-a*c*δθ² - a*c*δϕ² - a*s*δ²θ - b*c*δ²ϕ + b*s*δϕθ,
		)
	elseif rot_seq == :YZX
		return DCM(
			-a*c*δθ² - a*c*δϕ² - a*s*δ²θ - b*c*δ²ϕ + b*s*δϕθ,
			-a*c*d*δ²ϕ + a*c*e*δϕγ + a*d*s*δϕθ + b*c*d*δϕ² + δ²γ*(b*c*e + d*s) + δ²θ*(b*d*s + c*e) + δγ²*(b*c*d - e*s) + δθ²*(b*c*d - e*s) + δθγ*(-b*e*s + c*d),
			a*c*d*δϕγ + a*c*e*δ²ϕ - a*e*s*δϕθ - b*c*e*δϕ² + δ²γ*(b*c*d - e*s) + δ²θ*(-b*e*s + c*d) + δγ²*(-b*c*e - d*s) + δθ²*(-b*c*e - d*s) + δθγ*(-b*d*s - c*e),
			a*δ²ϕ - b*δϕ²,
			-a*d*δγ² - a*d*δϕ² - a*e*δ²γ - b*d*δ²ϕ + b*e*δϕγ,
			-a*d*δ²γ + a*e*δγ² + a*e*δϕ² + b*d*δϕγ + b*e*δ²ϕ,
			-a*c*δ²θ + a*s*δθ² + a*s*δϕ² + b*c*δϕθ + b*s*δ²ϕ,
			a*c*d*δϕθ + a*d*s*δ²ϕ - a*e*s*δϕγ - b*d*s*δϕ² + δ²γ*(-b*e*s + c*d) + δ²θ*(b*c*d - e*s) + δγ²*(-b*d*s - c*e) + δθ²*(-b*d*s - c*e) + δθγ*(-b*c*e - d*s),
			-a*c*e*δϕθ - a*d*s*δϕγ - a*e*s*δ²ϕ + b*e*s*δϕ² + δ²γ*(-b*d*s - c*e) + δ²θ*(-b*c*e - d*s) + δγ²*(b*e*s - c*d) + δθ²*(b*e*s - c*d) + δθγ*(-b*c*d + e*s),
		)
	elseif rot_seq == :YZY
		return DCM(
			-a*c*d*δϕ² - b*c*d*δ²ϕ + b*c*e*δϕγ + b*d*s*δϕθ + δ²γ*(-a*c*e - d*s) + δ²θ*(-a*d*s - c*e) + δγ²*(-a*c*d + e*s) + δθ²*(-a*c*d + e*s) + δθγ*(a*e*s - c*d),
			-a*c*δ²ϕ + a*s*δϕθ + b*c*δθ² + b*c*δϕ² + b*s*δ²θ,
			-a*c*e*δϕ² - b*c*d*δϕγ - b*c*e*δ²ϕ + b*e*s*δϕθ + δ²γ*(a*c*d - e*s) + δ²θ*(-a*e*s + c*d) + δγ²*(-a*c*e - d*s) + δθ²*(-a*c*e - d*s) + δθγ*(-a*d*s - c*e),
			a*d*δ²ϕ - a*e*δϕγ - b*d*δγ² - b*d*δϕ² - b*e*δ²γ,
			-a*δϕ² - b*δ²ϕ,
			a*d*δϕγ + a*e*δ²ϕ + b*d*δ²γ - b*e*δγ² - b*e*δϕ²,
			a*d*s*δϕ² + b*c*d*δϕθ + b*d*s*δ²ϕ - b*e*s*δϕγ + δ²γ*(a*e*s - c*d) + δ²θ*(-a*c*d + e*s) + δγ²*(a*d*s + c*e) + δθ²*(a*d*s + c*e) + δθγ*(a*c*e + d*s),
			a*c*δϕθ + a*s*δ²ϕ + b*c*δ²θ - b*s*δθ² - b*s*δϕ²,
			a*e*s*δϕ² + b*c*e*δϕθ + b*d*s*δϕγ + b*e*s*δ²ϕ + δ²γ*(-a*d*s - c*e) + δ²θ*(-a*c*e - d*s) + δγ²*(a*e*s - c*d) + δθ²*(a*e*s - c*d) + δθγ*(-a*c*d + e*s),
		)
	elseif rot_seq == :ZXY
		return DCM(
			-a*c*e*δϕθ - a*d*s*δϕγ - a*e*s*δ²ϕ + b*e*s*δϕ² + δ²γ*(-b*d*s - c*e) + δ²θ*(-b*c*e - d*s) + δγ²*(b*e*s - c*d) + δθ²*(b*e*s - c*d) + δθγ*(-b*c*d + e*s),
			-a*c*δ²θ + a*s*δθ² + a*s*δϕ² + b*c*δϕθ + b*s*δ²ϕ,
			a*c*d*δϕθ + a*d*s*δ²ϕ - a*e*s*δϕγ - b*d*s*δϕ² + δ²γ*(-b*e*s + c*d) + δ²θ*(b*c*d - e*s) + δγ²*(-b*d*s - c*e) + δθ²*(-b*d*s - c*e) + δθγ*(-b*c*e - d*s),
			a*c*d*δϕγ + a*c*e*δ²ϕ - a*e*s*δϕθ - b*c*e*δϕ² + δ²γ*(b*c*d - e*s) + δ²θ*(-b*e*s + c*d) + δγ²*(-b*c*e - d*s) + δθ²*(-b*c*e - d*s) + δθγ*(-b*d*s - c*e),
			-a*c*δθ² - a*c*δϕ² - a*s*δ²θ - b*c*δ²ϕ + b*s*δϕθ,
			-a*c*d*δ²ϕ + a*c*e*δϕγ + a*d*s*δϕθ + b*c*d*δϕ² + δ²γ*(b*c*e + d*s) + δ²θ*(b*d*s + c*e) + δγ²*(b*c*d - e*s) + δθ²*(b*c*d - e*s) + δθγ*(-b*e*s + c*d),
			-a*d*δ²γ + a*e*δγ² + a*e*δϕ² + b*d*δϕγ + b*e*δ²ϕ,
			a*δ²ϕ - b*δϕ²,
			-a*d*δγ² - a*d*δϕ² - a*e*δ²γ - b*d*δ²ϕ + b*e*δϕγ,
		)
	elseif rot_seq == :ZXZ
		return DCM(
			a*e*s*δϕ² + b*c*e*δϕθ + b*d*s*δϕγ + b*e*s*δ²ϕ + δ²γ*(-a*d*s - c*e) + δ²θ*(-a*c*e - d*s) + δγ²*(a*e*s - c*d) + δθ²*(a*e*s - c*d) + δθγ*(-a*c*d + e*s),
			a*d*s*δϕ² + b*c*d*δϕθ + b*d*s*δ²ϕ - b*e*s*δϕγ + δ²γ*(a*e*s - c*d) + δ²θ*(-a*c*d + e*s) + δγ²*(a*d*s + c*e) + δθ²*(a*d*s + c*e) + δθγ*(a*c*e + d*s),
			a*c*δϕθ + a*s*δ²ϕ + b*c*δ²θ - b*s*δθ² - b*s*δϕ²,
			-a*c*e*δϕ² - b*c*d*δϕγ - b*c*e*δ²ϕ + b*e*s*δϕθ + δ²γ*(a*c*d - e*s) + δ²θ*(-a*e*s + c*d) + δγ²*(-a*c*e - d*s) + δθ²*(-a*c*e - d*s) + δθγ*(-a*d*s - c*e),
			-a*c*d*δϕ² - b*c*d*δ²ϕ + b*c*e*δϕγ + b*d*s*δϕθ + δ²γ*(-a*c*e - d*s) + δ²θ*(-a*d*s - c*e) + δγ²*(-a*c*d + e*s) + δθ²*(-a*c*d + e*s) + δθγ*(a*e*s - c*d),
			-a*c*δ²ϕ + a*s*δϕθ + b*c*δθ² + b*c*δϕ² + b*s*δ²θ,
			a*d*δϕγ + a*e*δ²ϕ + b*d*δ²γ - b*e*δγ² - b*e*δϕ²,
			a*d*δ²ϕ - a*e*δϕγ - b*d*δγ² - b*d*δϕ² - b*e*δ²γ,
			-a*δϕ² - b*δ²ϕ,
		)
	elseif rot_seq == :ZYZ
		return DCM(
			-a*c*d*δϕ² - b*c*d*δ²ϕ + b*c*e*δϕγ + b*d*s*δϕθ + δ²γ*(-a*c*e - d*s) + δ²θ*(-a*d*s - c*e) + δγ²*(-a*c*d + e*s) + δθ²*(-a*c*d + e*s) + δθγ*(a*e*s - c*d),
			a*c*e*δϕ² + b*c*d*δϕγ + b*c*e*δ²ϕ - b*e*s*δϕθ + δ²γ*(-a*c*d + e*s) + δ²θ*(a*e*s - c*d) + δγ²*(a*c*e + d*s) + δθ²*(a*c*e + d*s) + δθγ*(a*d*s + c*e),
			a*c*δ²ϕ - a*s*δϕθ - b*c*δθ² - b*c*δϕ² - b*s*δ²θ,
			-a*d*s*δϕ² - b*c*d*δϕθ - b*d*s*δ²ϕ + b*e*s*δϕγ + δ²γ*(-a*e*s + c*d) + δ²θ*(a*c*d - e*s) + δγ²*(-a*d*s - c*e) + δθ²*(-a*d*s - c*e) + δθγ*(-a*c*e - d*s),
			a*e*s*δϕ² + b*c*e*δϕθ + b*d*s*δϕγ + b*e*s*δ²ϕ + δ²γ*(-a*d*s - c*e) + δ²θ*(-a*c*e - d*s) + δγ²*(a*e*s - c*d) + δθ²*(a*e*s - c*d) + δθγ*(-a*c*d + e*s),
			a*c*δϕθ + a*s*δ²ϕ + b*c*δ²θ - b*s*δθ² - b*s*δϕ²,
			-a*d*δ²ϕ + a*e*δϕγ + b*d*δγ² + b*d*δϕ² + b*e*δ²γ,
			a*d*δϕγ + a*e*δ²ϕ + b*d*δ²γ - b*e*δγ² - b*e*δϕ²,
			-a*δϕ² - b*δ²ϕ,
		)
	else
		throw(ArgumentError("The rotation sequence :$rot_seq is not valid."))
	end


end



function angle_to_δ³dcm(θ, rot_seq::Symbol)

    s, c = sincos(θ[1])

    a, b = (s, c).*θ[2].*θ[3].*3
    d, e = (s, c).*θ[2]^3
    f, g = (s, c).*θ[4]

    if rot_seq == :X 
        return DCM(0, 0, 0, 0, -f-b+d, -g+a+e, 0, g-a-e, -f-b+d) 

    elseif rot_seq == :Y 
        return DCM(-f-b+d, 0, g-a-e, 0, 0, 0, -g+a+e, 0, -f-b+d) 

    elseif rot_seq == :Z 
        return DCM(-f-b+d, -g+a+e, 0, g-a-e, -f-b+d, 0, 0, 0, 0) 

    else 
        throw(ArgumentError("rot_seq must be :X, :Y or :Z"))
    end

end