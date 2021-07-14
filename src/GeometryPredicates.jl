module GeometryPredicates
export two_sum, two_diff, two_prod,
       grow_expansion, expansion_sum, expansion_diff, expansion_scale

@inline function two_sum_tail(a::T, b::T, x::T) where {T}
    v = x - a
    return (a - (x - v)) + (b - v)
end

@inline function fast_two_sum_tail(a::T, b::T, x::T) where {T}
    v = x - a
    return b - v
end

@inline function two_sum(a::T, b::T) where {T}
    x = a + b
    return two_sum_tail(a, b, x), x
end

@inline function fast_two_sum(a::T, b::T) where {T}
    x = a + b
    return fast_two_sum_tail(a, b, x), x
end

@inline function two_diff_tail(a::T, b::T, x::T) where {T}
    b1 = a - x
    a1 = x + b1
    return (a - a1) + (b1 - b)
end

@inline function two_diff(a::T, b::T) where {T}
    x = a - b
    return two_diff_tail(a, b, x), x
end

@inline two_prod_tail(a::T, b::T, x::T) where {T} = fma(a, b, -x)

@inline function two_prod(a::T, b::T) where {T}
    x = a * b
    return two_prod_tail(a, b, x), x
end    

@inline grow_expansion(a::T, b::T) where {T} = two_sum(a, b)
@inline grow_expansion(a) = a

@generated function grow_expansion(e::NTuple{N,T}, b::T) where {N,T}
    h = [Symbol(:h, i) for i in 1:N]
    body = [:(($(h[i]), b) = two_sum(b, e[$i])) for i in 1:N]
    return quote
        $(body...)
        return ($(h...), b)
    end
end

@generated function expansion_sum(e::NTuple{N,T}, f::NTuple{M,T}) where {N,M,T}
    h = [Symbol(:h, i) for i in 1:(N+M)]
    body = [
        :(($(h[i:i+N]...),) = grow_expansion(($(h[i:i+N-1]...),),f[$i])) for i in 1:M
    ]
    return quote
        ($(h[1:N]...),) = e
        $(body...)
        return ($(h...),)
    end
end

@inline function expansion_diff(e::NTuple{2,T}, b::T) where {T}
    x1, x2 = two_diff(e[1], b)
    x2, x3 = two_sum(e[2], x2)
    return x1, x2, x3
end

@inline function expansion_diff(e::NTuple{2,T}, f::NTuple{2,T}) where {T}
    x1, x2, x3 = expansion_diff((e[2], e[1]), f[1])
    x2, x3, x4 = expansion_diff((x3, x2), f[2])
    return x1, x2, x3, x4
end

@generated function expansion_scale(e::NTuple{N,T}, b::T) where {N,T}
    h = [Symbol(:h, i) for i in 1:(2N-1)]
    Q = [Symbol(:Q, i) for i in 2:2N]
    t = [Symbol(:t, i) for i in 2:N]
    s = [Symbol(:s, i) for i in 2:N]
    body = [
        quote
            ($(t[i-1]), $(s[i-1])) = two_prod(e[$i], b)
            ($(h[2i-2]), $(Q[2i-2])) = two_sum($(t[i-1]), $(Q[2i-3]))
            ($(h[2i-1]), $(Q[2i-1])) = fast_two_sum($(s[i-1]), $(Q[2i-2]))
        end for i in 2:N
    ]
    return quote
        (h1, Q2) = two_prod(e[1], b)
        $(body...)
        return ($(h...), $(Q[end]))
    end
end

end # module
