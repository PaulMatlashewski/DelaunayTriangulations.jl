module DelaunayTriangulations

include("GeometryPredicates.jl")
using .GeometryPredicates
export orient_fast, orient_exact, orient, incircle_fast, incircle_exact, incircle

include("triangulation.jl")
export Triangulation, add_triangle!, delete_triangle!, adjacent, adjacent_to_vertex, orient, incircle

end # module
