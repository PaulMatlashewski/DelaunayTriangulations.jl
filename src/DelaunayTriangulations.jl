module DelaunayTriangulations

include("GeometryPredicates.jl")
using .GeometryPredicates
export orient_fast, orient_exact, orient, orient_adapt

end # module
