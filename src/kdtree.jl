import Random: randperm

struct Node
    value::Int
    left::Int
    right::Int
end

struct KDTree{T}
    values::Vector{NTuple{2,T}}
    nodes::Vector{Node}
end

function KDTree(values)
    n = length(values)
    nodes = [Node(0, 0, 0) for _ in 1:n]
    node_id = 1
    nodes[node_id] = Node(1, 0, 0)
    tree = KDTree(values, nodes)
    next_node = node_id + 1
    for value_id in randperm(n - 1)
        next_node = insert!(tree, value_id + 1, node_id, next_node, 1)
    end
    return tree
end

value(tree, value_id, dim) = tree.values[value_id][dim]
node(tree, node_id) = tree.nodes[node_id]

function insert!(tree, value_id, node_id, next_node, dim)
    root = node(tree, node_id)
    if value(tree, value_id, dim) > value(tree, root.value, dim)
        if root.right == 0
            tree.nodes[next_node] = Node(value_id, 0, 0)
            tree.nodes[node_id] = Node(root.value, root.left, next_node)
            next_node += 1
        else
            next_node = insert!(tree, value_id, root.right, next_node, dim % 2 + 1)
        end
    else
        if root.left == 0
            tree.nodes[next_node] = Node(value_id, 0, 0)
            tree.nodes[node_id] = Node(root.value, next_node, root.right)
            next_node += 1
        else
            next_node = insert!(tree, value_id, root.left, next_node, dim % 2 + 1)
        end
    end
    return next_node
end

function sort_values(tree)
    n = length(tree.values)
    permutation = zeros(Int, n)
    sort_values!(permutation, tree, 1, 1)
    return permutation
end

function sort_values!(permutation, tree, index, node_id)
    root = node(tree, node_id)
    if root.left != 0
        index = sort_values!(permutation, tree, index, root.left)
    end
    permutation[index] = root.value
    index += 1
    if root.right != 0
        index = sort_values!(permutation, tree, index, root.right)
    end
    return index
end

sort_vertices(vertices) = sort_values(KDTree(vertices))
