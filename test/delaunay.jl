using DelaunayTriangulations

@testset "Delaunay" begin
    vertices = [(0.0, 0.0), (0.5, 0.0), (1.0, 0.0),
                (0.0, 0.5), (0.5, 0.5), (1.0, 0.5),
                (0.0, 1.0), (0.5, 1.0), (1.0, 1.0)]
    triangulation = Triangulation(vertices)
    add_triangle!(triangulation, (1, 2, 5))
    add_triangle!(triangulation, (1, 5, 4))
    add_triangle!(triangulation, (2, 3, 6))
    add_triangle!(triangulation, (2, 6, 5))
    add_triangle!(triangulation, (4, 5, 7))
    add_triangle!(triangulation, (5, 8, 7))
    add_triangle!(triangulation, (5, 6, 8))
    add_triangle!(triangulation, (8, 6, 9))
    # Ghost triangles
    add_triangle!(triangulation, (2, 1, -1))
    add_triangle!(triangulation, (1, 4, -1))
    add_triangle!(triangulation, (4, 7, -1))
    add_triangle!(triangulation, (7, 8, -1))
    add_triangle!(triangulation, (8, 9, -1))
    add_triangle!(triangulation, (9, 6, -1))
    add_triangle!(triangulation, (6, 3, -1))
    add_triangle!(triangulation, (3, 2, -1))
    @testset "Insert Vertex" begin
        @testset "Inside Triangulation" begin
            # New vertex has index 10
            insert_vertex!(triangulation, (0.25, 0.01), (1, 2, 5))
            # Triangle (1, 2, 5) and (1, 5, 4) are split into
            # (1, 10, 4), (1, 2, 10), (10, 2, 5) and (10, 5, 4)
            @test adjacent(triangulation, (1, 2)) == 10
            @test adjacent(triangulation, (2, 10)) == 1
            @test adjacent(triangulation, (10, 1)) == 2
            @test adjacent(triangulation, (10, 2)) == 5
            @test adjacent(triangulation, (2, 5)) == 10
            @test adjacent(triangulation, (5, 10)) == 2
            @test adjacent(triangulation, (10, 5)) == 4
            @test adjacent(triangulation, (5, 4)) == 10
            @test adjacent(triangulation, (4, 10)) == 5
            @test adjacent(triangulation, (10, 4)) == 1
            @test adjacent(triangulation, (4, 1)) == 10
            @test adjacent(triangulation, (1, 10)) == 4
            # All other triangles are unchanged
            @test adjacent(triangulation, (2, 3)) == 6
            @test adjacent(triangulation, (3, 6)) == 2
            @test adjacent(triangulation, (6, 2)) == 3
            @test adjacent(triangulation, (6, 5)) == 2
            @test adjacent(triangulation, (5, 2)) == 6
            @test adjacent(triangulation, (2, 6)) == 5
            @test adjacent(triangulation, (4, 5)) == 7
            @test adjacent(triangulation, (5, 7)) == 4
            @test adjacent(triangulation, (7, 4)) == 5
            @test adjacent(triangulation, (5, 6)) == 8
            @test adjacent(triangulation, (6, 8)) == 5
            @test adjacent(triangulation, (8, 5)) == 6
            @test adjacent(triangulation, (6, 9)) == 8
            @test adjacent(triangulation, (9, 8)) == 6
            @test adjacent(triangulation, (8, 6)) == 9
            @test adjacent(triangulation, (2, 1)) == -1
            @test adjacent(triangulation, (1, 4)) == -1
            @test adjacent(triangulation, (4, 7)) == -1
            @test adjacent(triangulation, (7, 8)) == -1
            @test adjacent(triangulation, (8, 9)) == -1
            @test adjacent(triangulation, (9, 6)) == -1
            @test adjacent(triangulation, (6, 3)) == -1
            @test adjacent(triangulation, (3, 2)) == -1
        end
        @testset "Outside Triangulation" begin
            # New vertex has index 11
            insert_vertex!(triangulation, (1.1, 0.6), (9, 6, -1))
            # Ghost triangles (9, 6, -1) and (6, 3, -1) are split into
            # (6, 11, 9), (6, 3, 11), (9, 11, -1) and (11, 3, -1)
            @test adjacent(triangulation, (9, 6)) == 11
            @test adjacent(triangulation, (6, 11)) == 9
            @test adjacent(triangulation, (11, 9)) == 6
            @test adjacent(triangulation, (6, 3)) == 11
            @test adjacent(triangulation, (3, 11)) == 6
            @test adjacent(triangulation, (11, 6)) == 3
            @test adjacent(triangulation, (9, 11)) == -1
            @test adjacent(triangulation, (11, -1)) == 9
            @test adjacent(triangulation, (-1, 9)) == 11
        end
        @testset "Triangulation" begin
            vertices = [(rand(), rand()) for _ in 1:100]
            triangulation = delaunay_triangulation(vertices)
            for (u, v) in keys(triangulation.edges)
                w = adjacent(triangulation, (u, v))
                if u != -1 && v != -1 && w != -1
                    for q in 1:length(triangulation.vertices)
                        if q != u && q != v && q != w
                            @test incircle(triangulation, u, v, w, q) < 0.0
                        end
                    end
                end
            end
        end
    end
end
