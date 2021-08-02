function delaunay_triangulation(vertices::Vector{NTuple{2,T}}) where {T}
    n = length(vertices)
    ordered_vertices = biased_random_order(vertices)
    triangulation = Triangulation(ordered_vertices)
    triangle, w = initialize_triangulation!(triangulation)
    if w > 3
        # Add vertices skipped for initial triangulation
        for q in 3:w
            triangle = locate_vertex(triangulation, triangle, q)
            insert_vertex!(triangulation, q, triangle)
            u, v = adjacent_to_vertex(triangulation, q)
            triangle = q, u ,v
        end
    end
    for q in w:n
        triangle = locate_vertex(triangulation, triangle, q)
        insert_vertex!(triangulation, q, triangle)
        u, v = adjacent_to_vertex(triangulation, q)
        triangle = q, u ,v
    end
    return triangulation
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
        return (w, v, u), w
    else
        # Initial triangle
        add_triangle!(triangulation, (u, v, w))
        # Ghost triangles
        add_triangle!(triangulation, (v, u, -1))
        add_triangle!(triangulation, (u, w, -1))
        add_triangle!(triangulation, (w, v, -1))
        return (u, v, w), w
    end
end

function initial_vertices(vertices)
    n = length(vertices)
    u_index, u = 1, vertices[1]
    v_index, v = 2, vertices[2]
    for w_index in 3:n
        w = vertices[w_index]
        direction = orient(u, v, w)
        if direction != 0
            reverse = direction < 0
            return u_index, v_index, w_index, reverse
        end
    end
    if collinear_input
        error("all input vertices are collinear")
        return (0, 0, 0, false)
    end
end

function locate_vertex(triangulation, start_triangle, q)
    u, v, w = start_triangle
    # If the start triangle is a ghost triangle,
    # begin by walking into the triangulation.
    if u == -1
        u, v, w = w, v, adjacent(triangulation, (w, v))
    elseif v == -1
        u, v, w = u, w, adjacent(triangulation, (u, w))
    elseif w == -1
        u, v, w = v, u, adjacent(triangulation, (v, u))
    end
    while true
        # If the walk visits a ghost triangle,
        # the ghost triangle must contain the query point.
        if u == -1 || v == -1 || w == -1
            return u, v, w
        end
        if orient(triangulation, q, u, v) < 0
            u, v, w = v, u, adjacent(triangulation, (v, u))
            continue
        end
        if orient(triangulation, q, v, w) < 0
            u, v, w = w, v, adjacent(triangulation, (w, v))
            continue
        end
        if orient(triangulation, q, w, u) < 0
            u, v, w = u, w, adjacent(triangulation, (u, w))
            continue
        end
        return u, v, w
    end
end

function insert_vertex!(triangulation::Triangulation{T}, vertex::NTuple{2,T}, triangle::NTuple{3,Int}) where {T}
    push!(triangulation.vertices, vertex)
    insert_vertex!(triangulation, length(triangulation.vertices), triangle)
end

function insert_vertex!(triangulation::Triangulation, u::Int, triangle::NTuple{3,Int})
    delete_triangle!(triangulation, triangle)
    v, w, x = triangle
    dig_cavity!(triangulation, u, (v, w))
    dig_cavity!(triangulation, u, (w, x))
    dig_cavity!(triangulation, u, (x, v))
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
end
