using DelaunayTriangulations:
    KDTree,
    Node,
    spatial_order,
    biased_random_order

function test_binary_search_property(tree, node_id, dim)
    root = tree.nodes[node_id]
    if root.left != 0
        left = tree.nodes[root.left]
        @test root.value[dim] >= left.value[dim]
        test_binary_search_property(tree, root.left, dim % 2 + 1)
    end
    if root.right != 0
        right = tree.nodes[root.right]
        @test root.value[dim] < right.value[dim]
        test_binary_search_property(tree, root.right, dim % 2 + 1)
    end
end

@testset "KDTree" begin
    # Check edge case of one vertex in the tree
    vertices = [(rand(), rand())]
    tree = KDTree(vertices)
    @test tree.nodes == [Node(vertices[1], 1, 0, 0)]
    # Check binary search property
    tree = KDTree([(rand(), rand()) for _ in 1:10])
    test_binary_search_property(tree, 1, 1)
    # Check sorting
    vertices = [(0.5185, 0.0347), (0.5216, 0.4219),
                (0.6828, 0.9237), (0.4302, 0.4148),
                (0.1816, 0.2890), (0.9312, 0.6561),
                (0.5897, 0.3000), (0.8778, 0.0046),
                (0.5463, 0.4520), (0.2716, 0.1995)]
    nodes = [Node(vertices[1], 1, 8, 2), Node(vertices[2], 2, 3, 4),
             Node(vertices[8], 8, 7, 0), Node(vertices[6], 6, 5, 0),
             Node(vertices[9], 9, 0, 6), Node(vertices[3], 3, 0, 0),
             Node(vertices[7], 7, 0, 0), Node(vertices[4], 4, 9, 0),
             Node(vertices[5], 5, 0, 10), Node(vertices[10], 10, 0, 0)]
    vertex_order = spatial_order(KDTree(nodes))
    @test vertex_order == [4, 7, 9, 3, 1, 10, 5, 6, 8, 2]
end
