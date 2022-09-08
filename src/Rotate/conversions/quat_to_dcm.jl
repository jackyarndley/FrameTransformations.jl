export quat_to_dcm

"""
    quat_to_dcm(q::Quaternion)

Convert the quaternion `q` to a Direction Cosine Matrix (DCM).
"""
function quat_to_dcm(q::Quaternion)
    # Auxiliary variables.
    q0 = q.q0
    q1 = q.q1
    q2 = q.q2
    q3 = q.q3

    return DCM(
        q0^2 + q1^2 - q2^2 - q3^2,   2(q1 * q2 + q0 * q3)   ,   2(q1 * q3 - q0 * q2),
          2(q1 * q2 - q0 * q3)   , q0^2 - q1^2 + q2^2 - q3^2,   2(q2 * q3 + q0 * q1),
          2(q1 * q3 + q0 * q2)   ,   2(q2 * q3 - q0 * q1)   , q0^2 - q1^2 - q2^2 + q3^2
    )'
end
