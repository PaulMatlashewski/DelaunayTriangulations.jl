module GeometryPredicates
export orient_fast, orient_exact, orient, incircle_fast, incircle_exact, incircle

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
const EPS_INCIRCLE_A_FLOAT64 = (10.0 + 96.0 * eps(Float64)) * eps(Float64)
const EPS_INCIRCLE_A_FLOAT32 = (10.0 + 96.0 * eps(Float32)) * eps(Float32)
const EPS_INCIRCLE_A_FLOAT16 = (10.0 + 96.0 * eps(Float16)) * eps(Float16)

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

function incircle_fast(a::NTuple{2,T}, b::NTuple{2,T}, c::NTuple{2,T}, d::NTuple{2,T}) where {T}
    ad1, ad2 = a[1] - d[1], a[2] - d[2]
    bd1, bd2 = b[1] - d[1], b[2] - d[2]
    cd1, cd2 = c[1] - d[1], c[2] - d[2]
    ab_det = ad1 * bd2 - bd1 * ad2
    bc_det = bd1 * cd2 - cd1 * bd2
    ca_det = cd1 * ad2 - ad1 * cd2
    a_lift = ad1 * ad1 + ad2 * ad2
    b_lift = bd1 * bd1 + bd2 * bd2
    c_lift = cd1 * cd1 + cd2 * cd2
    return a_lift * bc_det + b_lift * ca_det + c_lift * ab_det
end

function incircle_exact(a::NTuple{2,T}, b::NTuple{2,T}, c::NTuple{2,T}, d::NTuple{2,T}) where {T}
    ab_1, ab_2 = two_prod(a[1], b[2]), two_prod(b[1], a[2])
    bc_1, bc_2 = two_prod(b[1], c[2]), two_prod(c[1], b[2])
    cd_1, cd_2 = two_prod(c[1], d[2]), two_prod(d[1], c[2])
    da_1, da_2 = two_prod(d[1], a[2]), two_prod(a[1], d[2])
    ac_1, ac_2 = two_prod(a[1], c[2]), two_prod(c[1], a[2])
    bd_1, bd_2 = two_prod(b[1], d[2]), two_prod(d[1], b[2])

    ab = expansion_diff(ab_1, ab_2)
    bc = expansion_diff(bc_1, bc_2)
    cd = expansion_diff(cd_1, cd_2)
    da = expansion_diff(da_1, da_2)
    ac = expansion_diff(ac_1, ac_2)
    bd = expansion_diff(bd_1, bd_2)

    ca = (-ac[1], -ac[2], -ac[3], -ac[4])
    db = (-bd[1], -bd[2], -bd[3], -bd[4])

    cda = expansion_sum(expansion_sum(cd, da), ac)
    dab = expansion_sum(expansion_sum(da, ab), bd)
    abc = expansion_sum(expansion_sum(ab, bc), ca)
    bcd = expansion_sum(expansion_sum(bc, cd), db)

    adet_1 = expansion_scale(expansion_scale(bcd, a[1]), a[1])
    adet_2 = expansion_scale(expansion_scale(bcd, a[2]), a[2])
    bdet_1 = expansion_scale(expansion_scale(cda, b[1]), -b[1])
    bdet_2 = expansion_scale(expansion_scale(cda, b[2]), -b[2])
    cdet_1 = expansion_scale(expansion_scale(dab, c[1]), c[1])
    cdet_2 = expansion_scale(expansion_scale(dab, c[2]), c[2])
    ddet_1 = expansion_scale(expansion_scale(abc, d[1]), -d[1])
    ddet_2 = expansion_scale(expansion_scale(abc, d[2]), -d[2])

    adet = expansion_sum(adet_1, adet_2)
    bdet = expansion_sum(bdet_1, bdet_2)
    cdet = expansion_sum(cdet_1, cdet_2)
    ddet = expansion_sum(ddet_1, ddet_2)

    det = expansion_sum(expansion_sum(adet, bdet), expansion_sum(cdet, ddet))
    return det[end]
end

@inline incircle_test_A(A::Float64, alpha::Float64) = abs(A) >= EPS_INCIRCLE_A_FLOAT64 * alpha
@inline incircle_test_A(A::Float32, alpha::Float32) = abs(A) >= EPS_INCIRCLE_A_FLOAT32 * alpha
@inline incircle_test_A(A::Float16, alpha::Float16) = abs(A) >= EPS_INCIRCLE_A_FLOAT16 * alpha

# TODO: Add adaptive method instead of falling back to exact
function incircle(a::NTuple{2,T}, b::NTuple{2,T}, c::NTuple{2,T}, d::NTuple{2,T}) where {T}
    ad1, ad2 = a[1] - d[1], a[2] - d[2]
    bd1, bd2 = b[1] - d[1], b[2] - d[2]
    cd1, cd2 = c[1] - d[1], c[2] - d[2]
    adbd, bdad = ad1 * bd2, bd1 * ad2
    bdcd, cdbd = bd1 * cd2, cd1 * bd2
    cdad, adcd = cd1 * ad2, ad1 * cd2
    a_lift = ad1 * ad1 + ad2 * ad2
    b_lift = bd1 * bd1 + bd2 * bd2
    c_lift = cd1 * cd1 + cd2 * cd2
    A = (bdcd - cdbd) * a_lift + (cdad - adcd) * b_lift + (adbd - bdad) * c_lift
    alpha = (
        (abs(bdcd) + abs(cdbd)) * a_lift + 
        (abs(cdad) + abs(adcd)) * b_lift + 
        (abs(adbd) + abs(bdad)) * c_lift
    )
    if incircle_test_A(A, alpha)
        return A
    else
        return incircle_exact(a, b, c, d)
    end
end

end # module
