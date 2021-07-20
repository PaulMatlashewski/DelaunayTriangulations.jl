using DelaunayTriangulations: KDTree, Node, sort_vertices, sort_values, value, node

function test_binary_search_property(tree, node_id, dim)
    root = node(tree, node_id)
    if root.left != 0
        left = node(tree, root.left)
        @test value(tree, root.value, dim) >= value(tree, left.value, dim)
        test_binary_search_property(tree, root.left, dim % 2 + 1)
    end
    if root.right != 0
        right = node(tree, root.right)
        @test value(tree, root.value, dim) < value(tree, right.value, dim)
        test_binary_search_property(tree, root.right, dim % 2 + 1)
    end
end

@testset "KDTree" begin
    # Check edge case of one vertex in the tree
    vertices = [(rand(), rand())]
    tree = KDTree(vertices)
    @test tree.values == vertices
    @test tree.nodes == [Node(1, 0, 0)]
    # Check binary search property
    tree = KDTree([(rand(), rand()) for _ in 1:10])
    test_binary_search_property(tree, 1, 1)
    # Check permutation order
    vertices = [(0.5185, 0.0347), (0.5216, 0.4219),
                (0.6828, 0.9237), (0.4302, 0.4148),
                (0.1816, 0.2890), (0.9312, 0.6561),
                (0.5897, 0.3000), (0.8778, 0.0046),
                (0.5463, 0.4520), (0.2716, 0.1995)]
    nodes = [Node(1, 8, 2), Node(2, 3, 4), Node(8, 7, 0), Node(6, 5, 0), Node(9, 0, 6),
             Node(3, 0, 0), Node(7, 0, 0), Node(4, 9, 0), Node(5, 0, 10), Node(10, 0, 0)]
    tree = KDTree(vertices, nodes)
    permutation = sort_values(tree)
    @test permutation == [5, 10, 4, 1, 7, 8, 2, 9, 3, 6]
    @test length(sort_vertices(vertices)) == length(permutation)
end
