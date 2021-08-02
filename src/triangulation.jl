import .GeometryPredicates: orient, incircle

struct Triangulation{T}
    vertices::Vector{NTuple{2,T}}
    edges::Dict{NTuple{2,Int},Int}
end

function Triangulation(vertices::Vector{NTuple{2,T}}) where {T}
    return Triangulation(vertices, Dict{NTuple{2,Int},Int}())
end

function Triangulation(T::DataType, n::Int)
    return Triangulation(Vector{NTuple{2,T}}(undef, n), Dict{NTuple{2,Int},Int}())
end

function add_triangle!(triangulation::Triangulation, triangle::NTuple{3,Int})
    u, v, w = triangle
    edges = keys(triangulation.edges)
    if ((u, v) in edges) || ((v, w) in edges) || ((w, u) in edges)
        return triangle
    end
    triangulation.edges[(u, v)] = w
    triangulation.edges[(v, w)] = u
    triangulation.edges[(w, u)] = v
    return triangle
end

function delete_triangle!(triangulation::Triangulation, triangle::NTuple{3,Int})
    u, v, w = triangle
    edges = keys(triangulation.edges)
    if ((u, v) in edges) && ((v, w) in edges) && ((w, u) in edges)
        delete!(triangulation.edges, (u, v))
        delete!(triangulation.edges, (v, w))
        delete!(triangulation.edges, (w, u))
    end
    return
end

function adjacent(triangulation::Triangulation, edge::NTuple{2,Int})
    return get(triangulation.edges, edge, 0)
end

function adjacent_to_vertex(triangulation::Triangulation, u::Int)
    for (key, value) in triangulation.edges
        if value == u
            return key
        end
    end
    return (0, 0)
end

function GeometryPredicates.orient(triangulation::Triangulation, u::Int, v::Int, w::Int)
    return orient(triangulation.vertices[u], triangulation.vertices[v],
                  triangulation.vertices[w])
end

function GeometryPredicates.incircle(triangulation::Triangulation, u::Int, v::Int, w::Int, q::Int)
    return incircle(triangulation.vertices[u], triangulation.vertices[v],
                    triangulation.vertices[w], triangulation.vertices[q])
end

function GeometryPredicates.incircle(triangulation::Triangulation, u::Int, v::Int, w::Int)
    return incircle(triangulation.vertices[u], triangulation.vertices[v],
                    triangulation.vertices[w])
end
