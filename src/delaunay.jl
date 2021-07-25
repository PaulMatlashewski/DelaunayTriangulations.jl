function delaunay_triangulation(vertices::Vector{NTuple{2,T}}) where {T}
    n = length(vertices)
    ordered_vertices = biased_random_order(vertices)
    triangulation = Triangulation(ordered_vertices)
    w = initialize_triangulation!(triangulation)
    if w > 3
        # Add vertices skipped for initial triangulation
        for i in 3:w
            triangle = locate_vertex(triangulation, i)
            insert_vertex!(triangulation, i, triangle)
        end
    end
    for i in w:n
        triangle = locate_vertex(triangulation, i)
        insert_vertex!(triangulation, i, triangle)
    end
    return triangulation
end

function initial_vertices(vertices)
    u_index, u = 1, vertices[1]
    v_index, v = 2, vertices[2]
    w_index, w = find_independent_vertex(vertices, u, v, 3:n)
    reverse = false
    if orient(u, v, w) < 0
        reverse = true
    end
    return u_index, v_index, w_index, reverse
end

function initialize_triangulation!(triangulation)
    u, v, w, reverse = initial_vertices(triangulation.vertices)
    if reverse
        # Initial triangle
        add_triangle!(triangulation, (w, v, u))
        # Ghost triangles
        add_triangle!(triangulation, (u, v, -1))
        add_triangle!(triangulation, (v, w, -1))
        add_triangle!(triangulation, (w, u, -1))
    else
        # Initial triangle
        add_triangle!(triangulation, (u, v, w))
        # Ghost triangles
        add_triangle!(triangulation, (v, u, -1))
        add_triangle!(triangulation, (u, w, -1))
        add_triangle!(triangulation, (w, v, -1))
    end
    return w
end

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
