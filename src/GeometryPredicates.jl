module GeometryPredicates
export orient_fast, orient_exact, orient, orient_adapt

include("adaptive_arithmetic.jl")

const EPS_ORIENT_A_FLOAT64 = (3.0 + 16.0 * eps(Float64)) * eps(Float64)
const EPS_ORIENT_B_FLOAT64 = (2.0 + 12.0 * eps(Float64)) * eps(Float64)
const EPS_ORIENT_C1_FLOAT64 = (3.0 + 8.0 * eps(Float64)) * eps(Float64)
const EPS_ORIENT_C2_FLOAT64 = (9.0 + 64.0 * eps(Float64)) * eps(Float64) * eps(Float64)
const EPS_ORIENT_A_FLOAT32 = (3.0 + 16.0 * eps(Float32)) * eps(Float32)
const EPS_ORIENT_B_FLOAT32 = (2.0 + 12.0 * eps(Float32)) * eps(Float32)
const EPS_ORIENT_C1_FLOAT32 = (3.0 + 8.0 * eps(Float32)) * eps(Float32)
const EPS_ORIENT_C2_FLOAT32 = (9.0 + 64.0 * eps(Float32)) * eps(Float32) * eps(Float32)
const EPS_ORIENT_A_FLOAT16 = (3.0 + 16.0 * eps(Float16)) * eps(Float16)
const EPS_ORIENT_B_FLOAT16 = (2.0 + 12.0 * eps(Float16)) * eps(Float16)
const EPS_ORIENT_C1_FLOAT16 = (3.0 + 8.0 * eps(Float16)) * eps(Float16)
const EPS_ORIENT_C2_FLOAT16 = (9.0 + 64.0 * eps(Float16)) * eps(Float16) * eps(Float16)

function orient_fast(a::NTuple{2,T}, b::NTuple{2,T}, c::NTuple{2,T}) where {T}
    ac1, ac2 = a[1] - c[1], a[2] - c[2]
    bc1, bc2 = b[1] - c[1], b[2] - c[2]
    return ac1 * bc2 - ac2 * bc1
end

function orient_exact(a::NTuple{2,T}, b::NTuple{2,T}, c::NTuple{2,T}) where {T}
    ab = two_prod(a[1], b[2])
    ac = two_prod(a[1], c[2])
    bc = two_prod(b[1], c[2])
    ba = two_prod(b[1], a[2])
    ca = two_prod(c[1], a[2])
    cb = two_prod(c[1], b[2])
    u = expansion_sum(expansion_diff(ab, ac), expansion_diff(bc, ba))
    v = expansion_sum(u, expansion_diff(ca, cb))
    return v[end]
end

@inline orient_test_A(A::Float64, detsum::Float64) = abs(A) >= EPS_ORIENT_A_FLOAT64 * detsum
@inline orient_test_A(A::Float32, detsum::Float32) = abs(A) >= EPS_ORIENT_A_FLOAT32 * detsum
@inline orient_test_A(A::Float16, detsum::Float16) = abs(A) >= EPS_ORIENT_A_FLOAT16 * detsum
@inline orient_test_B(B::Float64, detsum::Float64) = abs(B) >= EPS_ORIENT_B_FLOAT64 * detsum
@inline orient_test_B(B::Float32, detsum::Float32) = abs(B) >= EPS_ORIENT_B_FLOAT32 * detsum
@inline orient_test_B(B::Float16, detsum::Float16) = abs(B) >= EPS_ORIENT_B_FLOAT16 * detsum
@inline function orient_test_C(C::Float64, B::Float64, detsum::Float64)
    return abs(C) >= EPS_ORIENT_C1_FLOAT64 * abs(B) + EPS_ORIENT_C2_FLOAT64 * detsum
end
@inline function orient_test_C(C::Float32, B::Float32, detsum::Float32)
    return abs(C) >= EPS_ORIENT_C1_FLOAT32 * abs(B) + EPS_ORIENT_C2_FLOAT32 * detsum
end
@inline function orient_test_C(C::Float16, B::Float16, detsum::Float16)
    return abs(C) >= EPS_ORIENT_C1_FLOAT16 * abs(B) + EPS_ORIENT_C2_FLOAT16 * detsum
end

function orient(a::NTuple{2,T}, b::NTuple{2,T}, c::NTuple{2,T}) where {T}
    ac1, ac2 = a[1] - c[1], a[2] - c[2]
    bc1, bc2 = b[1] - c[1], b[2] - c[2]
    det1 = ac1 * bc2
    det2 = ac2 * bc1
    det = det1 - det2
    detsum = det1 + det2
    # The result is only questionable when det is O(epsilon)
    if det1 > 0.0
        if det2 <= 0.0
            return det
        end
    elseif det1 < 0.0
        if det2 >= 0.01
            return det
        else
            detsum *= -1.0
        end
    else
        return det
    end
    if orient_test_A(det, detsum)
        return det
    else
        return orient_adapt(a, b, c, ac1, ac2, bc1, bc2, detsum)
    end
end

function orient_adapt(a::NTuple{2,T}, b::NTuple{2,T}, c::NTuple{2,T},
                      ac1::T, ac2::T, bc1::T, bc2::T, detsum::T) where {T}
    det1 = two_prod(ac1, bc2)
    det2 = two_prod(ac2, bc1)
    B = expansion_diff(det1, det2)
    B_prime = approximate(B)
    if orient_test_B(B_prime, detsum)
        return B_prime
    end

    ac1_error = two_diff_tail(a[1], c[1], ac1)
    ac2_error = two_diff_tail(a[2], c[2], ac2)
    bc1_error = two_diff_tail(b[1], c[1], bc1)
    bc2_error = two_diff_tail(b[2], c[2], bc2)
    if (ac1_error == 0.0) && (ac2_error == 0.0) && (bc1_error == 0.0) && (bc2_error == 0.0)
        return B_prime
    end

    C = B_prime + (ac1 * bc2_error + bc2 * ac1_error) - (ac2 * bc1_error + bc1 * ac2_error)
    if orient_test_C(C, B_prime, detsum)
        return C
    end

    s = two_prod(ac1_error, bc2)
    t = two_prod(ac2_error, bc1)
    u = expansion_diff(s, t)
    C1 = expansion_sum(B, u)

    s = two_prod(ac1, bc2_error)
    t = two_prod(ac2, bc1_error)
    u = expansion_diff(s, t)
    C2 = expansion_sum(C1, u)

    s = two_prod(ac1_error, bc2_error)
    t = two_prod(ac2_error, bc1_error)
    u = expansion_diff(s, t)
    D = expansion_sum(C2, u)

    return D[end]
end

end # module
