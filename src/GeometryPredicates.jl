module GeometryPredicates
export two_sum, grow_expansion, expansion_sum

@inline function two_sum_tail(a::T, b::T, x::T) where {T}
    v = x - a
    return (a - (x - v)) + (b - v)
end

@inline function two_sum(a::T, b::T) where {T}
    x = a + b
    return two_sum_tail(a, b, x), x
end

grow_expansion(a::T, b::T) where {T} = two_sum(a, b)
grow_expansion(a) = a

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

end # module
