import Random: shuffle

struct Node{T}
    value::NTuple{2,T}
    value_id::Int
    left::Int
    right::Int
end

struct KDTree{T}
    nodes::Vector{Node{T}}
end

function KDTree(values::Vector{NTuple{2,T}}) where {T}
    n = length(values)
    nodes = [Node((0.0, 0.0), 0, 0, 0) for _ in 1:n]
    node_id = 1
    nodes[node_id] = Node(values[node_id], node_id, 0, 0)
    tree = KDTree(nodes)
    next_node = node_id + 1
    for value_id in 2:n
        next_node = insert!(tree, values[value_id], value_id, node_id, next_node, 1)
    end
    return tree
end

function insert!(tree, value, value_id, node_id, next_node, dim)
    root = tree.nodes[node_id]
    if value[dim] > root.value[dim]
        if root.right == 0
            tree.nodes[next_node] = Node(value, value_id, 0, 0)
            tree.nodes[node_id] = Node(root.value, root.value_id, root.left, next_node)
            next_node += 1
        else
            next_node = insert!(tree, value, value_id, root.right, next_node, dim % 2 + 1)
        end
    else
        if root.left == 0
            tree.nodes[next_node] = Node(value, value_id, 0, 0)
            tree.nodes[node_id] = Node(root.value, root.value_id, next_node, root.right)
            next_node += 1
        else
            next_node = insert!(tree, value, value_id, root.left, next_node, dim % 2 + 1)
        end
    end
    return next_node
end

function spatial_order(tree)
    n = length(tree.nodes)
    order = zeros(Int, n)
    spatial_order!(order, tree, 1, 1)
    return order
end

function spatial_order!(order, tree, index, node_id)
    root = tree.nodes[node_id]
    if root.left != 0
        index = spatial_order!(order, tree, index, root.left)
    end
    order[root.value_id] = index
    index += 1
    if root.right != 0
        index = spatial_order!(order, tree, index, root.right)
    end
    return index
end

function biased_random_order(vertices)
    n = length(vertices)
    m = Int(ceil(log(2, n))) # Number of insertion rounds
    vertices = shuffle(vertices)
    order = spatial_order(KDTree(vertices))
    reverse = true
    start = 1
    for i in 1:m - 1
        count = Int(max(1, floor(n/2^(m + 1 - i))))
        stop = start + count - 1
        order_round = @view order[start:stop]
        perm = sortperm(order_round; rev=reverse)
        vertices[start:stop] .= @view vertices[start:stop][perm]
        reverse = !reverse
        start = stop + 1
    end
    order_round = @view order[start:n]
    perm = sortperm(order_round; rev=reverse)
    vertices[start:n] .= @view vertices[start:n][perm]
    return vertices
end
