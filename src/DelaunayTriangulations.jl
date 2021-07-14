module DelaunayTriangulations

include("GeometryPredicates.jl")
using .GeometryPredicates
export orient_fast, orient_exact, orient, incircle_fast, incircle_exact, incircle

end # module
