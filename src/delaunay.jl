function insert_vertex!(triangulation::Triangulation{T}, vertex::NTuple{2,T}, triangle::NTuple{3,Int}) where {T}
    push!(triangulation.vertices, vertex)
    insert_vertex!(triangulation, length(triangulation.vertices), triangle)
    return
end

function insert_vertex!(triangulation::Triangulation, u::Int, triangle::NTuple{3,Int})
    delete_triangle!(triangulation, triangle)
    v, w, x = triangle
    dig_cavity!(triangulation, u, (v, w))
    dig_cavity!(triangulation, u, (w, x))
    dig_cavity!(triangulation, u, (x, v))
    return
end

function dig_cavity!(triangulation::Triangulation{T}, u, edge) where {T}
    v, w = edge
    x = adjacent(triangulation, (w, v))
    if x == 0 # Triangle has already been deleted
        return
    end
    if x == -1
    end
    if x == -1
        u_incircle = incircle(triangulation, u, w, v)
    elseif v == -1
        u_incircle = incircle(triangulation, u, x, w)
    elseif w == -1
        u_incircle = incircle(triangulation, u, v, x)
    else
        u_incircle = incircle(triangulation, u, v, w, x)
    end
    if u_incircle > zero(T)
        # uvw and uvx are not Delaunay
        delete_triangle!(triangulation, (w, v, x))
        dig_cavity!(triangulation, u, (v, x))
        dig_cavity!(triangulation, u, (x, w))
        return
    end
    add_triangle!(triangulation, (u, v, w))
    return
end
