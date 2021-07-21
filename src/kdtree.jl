import Random: shuffle

struct Node{T}
    value::NTuple{2,T}
    left::Int
    right::Int
end

struct KDTree{T}
    nodes::Vector{Node{T}}
end

function KDTree(values::Vector{NTuple{2,T}}) where {T}
    values = shuffle(values)
    n = length(values)
    nodes = [Node((0.0, 0.0), 0, 0) for _ in 1:n]
    node_id = 1
    nodes[node_id] = Node(values[node_id], 0, 0)
    tree = KDTree(nodes)
    next_node = node_id + 1
    for value_id in 2:n
        next_node = insert!(tree, values[value_id], node_id, next_node, 1)
    end
    return tree
end

function insert!(tree, value, node_id, next_node, dim)
    root = tree.nodes[node_id]
    if value[dim] > root.value[dim]
        if root.right == 0
            tree.nodes[next_node] = Node(value, 0, 0)
            tree.nodes[node_id] = Node(root.value, root.left, next_node)
            next_node += 1
        else
            next_node = insert!(tree, value, root.right, next_node, dim % 2 + 1)
        end
    else
        if root.left == 0
            tree.nodes[next_node] = Node(value, 0, 0)
            tree.nodes[node_id] = Node(root.value, next_node, root.right)
            next_node += 1
        else
            next_node = insert!(tree, value, root.left, next_node, dim % 2 + 1)
        end
    end
    return next_node
end

function sort_vertices(vertices)
    sorted_vertices = copy(vertices)
    tree = KDTree(vertices)
    sort_values!(sorted_vertices, tree, 1, 1)
    return sorted_vertices
end

function sort_vertices!(vertices)
    tree = KDTree(vertices)
    sort_values!(vertices, tree, 1, 1)
end

function sort_values!(values, tree, index, node_id)
    root = tree.nodes[node_id]
    if root.left != 0
        index = sort_values!(values, tree, index, root.left)
    end
    values[index] = root.value
    index += 1
    if root.right != 0
        index = sort_values!(values, tree, index, root.right)
    end
    return index
end
