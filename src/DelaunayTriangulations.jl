module DelaunayTriangulations

include("GeometryPredicates.jl")
using .GeometryPredicates
export two_sum, two_sum_tail, two_diff, two_diff_tail, two_prod, two_prod_tail,
       expansion_sum, expansion_diff, expansion_scale, approximate

end # module
