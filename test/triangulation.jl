using DelaunayTriangulations

@testset "Triangulation" begin
    vertices = [(0.0, 0.0), (0.5, 0.0), (1.0, 0.0),
                (0.0, 0.5), (0.5, 0.5), (1.0, 0.5),
                (0.0, 1.0), (0.5, 1.0), (1.0, 1.0)]
    triangulation = Triangulation(vertices)
    @testset "Triangulation Construction" begin
        add_triangle!(triangulation, (1, 2, 5))
        @test ((1, 2) => 5) in triangulation.edges
        @test ((2, 5) => 1) in triangulation.edges
        @test ((5, 1) => 2) in triangulation.edges
        @test length(triangulation.edges) == 3
        @test adjacent(triangulation, (1, 2)) == 5
        @test adjacent(triangulation, (2, 5)) == 1
        @test adjacent(triangulation, (5, 1)) == 2
        @test adjacent(triangulation, (1, 5)) == 0 # No triangle containing edge
    end
    @testset "Add Triangles" begin
        # If (u,v,w) is an existing triangle, triangles (u,v,x),
        # (x,v,w), and (u,x,w) are rejected for any x
        add_triangle!(triangulation, (1, 2, 6))
        add_triangle!(triangulation, (3, 2, 5))
        add_triangle!(triangulation, (1, 4, 5))
        @test ((1, 2) => 5) in triangulation.edges
        @test ((2, 5) => 1) in triangulation.edges
        @test ((5, 1) => 2) in triangulation.edges
        @test length(triangulation.edges) == 3
        add_triangle!(triangulation, (1, 5, 4))
        add_triangle!(triangulation, (2, 3, 6))
        add_triangle!(triangulation, (2, 6, 5))
        add_triangle!(triangulation, (4, 5, 7))
        add_triangle!(triangulation, (5, 8, 7))
        add_triangle!(triangulation, (5, 6, 8))
        add_triangle!(triangulation, (8, 6, 9))
        @test length(triangulation.edges) == 24
        @test adjacent_to_vertex(triangulation, 1) in Set([(2, 5), (5, 4)])
        @test adjacent_to_vertex(triangulation, 2) in Set([(5, 1), (6, 5)])
        @test adjacent_to_vertex(triangulation, 3) == (6, 2)
        @test adjacent_to_vertex(triangulation, 9) == (8, 6)
    end
    @testset "Geometry Predicates" begin
        push!(triangulation.vertices, (0.25, 0.0001))
        push!(triangulation.vertices, (0.25, 0.0))
        @testset "incircle" begin
            @test incircle(triangulation, 1, 2, 5, 10) > 0.0
            @test incircle(triangulation, 10, 1, 2) > 0.0
            @test incircle(triangulation, 10, 2, 1) < 0.0
            @test incircle(triangulation, 3, 1, 2) == 0.0
            @test incircle(triangulation, 11, 1, 2) == 1.0
        end
        @testset "orient" begin
            @test orient(triangulation, 3, 1, 2) == 0.0
            @test orient(triangulation, 5, 1, 2) > 0.0
            @test orient(triangulation, 4, 5, 2) < 0.0
        end
    end
    @testset "DeleteTriangles" begin
        delete_triangle!(triangulation, (5, 2, 1)) # Triangle not in triangulation, nothing should change
        @test length(triangulation.edges) == 24
        @test adjacent(triangulation, (5, 1)) == 2
        delete_triangle!(triangulation, (2, 5, 1))
        @test length(triangulation.edges) == 21
        @test adjacent(triangulation, (5, 1)) == 0
        delete_triangle!(triangulation, (5, 2, 6))
        @test adjacent_to_vertex(triangulation, 2) == (3, 6)
        delete_triangle!(triangulation, (3, 6, 2))
        @test adjacent_to_vertex(triangulation, 2) == (0, 0)
        delete_triangle!(triangulation, (5, 4, 1))
        delete_triangle!(triangulation, (9, 8, 6))
        delete_triangle!(triangulation, (8, 5, 6))
        delete_triangle!(triangulation, (4, 5, 7))
        delete_triangle!(triangulation, (8, 7, 5))
        @test length(triangulation.edges) == 0
    end
end
