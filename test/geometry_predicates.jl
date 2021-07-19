using DelaunayTriangulations

@testset "Geometry Predicates" begin
    @testset "Orient" begin
        @testset "Orient Fast" begin
            @test orient_fast((0.1, 0.1), (0.2, 0.2), (0.3, 0.3)) == 0.0
            @test orient_fast((0.0, 3.14159e23), (0.0, 42.4242e11), (0.0, 1.234567e2)) == 0.0
            @test orient_fast((-1.0e-7, -1.0e-7), (2.0, 2.0), (3.0e7, 3.0e7)) == 0.0
            @test orient_fast((Float32(0.1), Float32(0.1)), (Float32(0.2), Float32(0.2)), (Float32(0.3), Float32(0.3))) == Float32(0.0)
            @test orient_fast((Float16(0.1), Float16(0.1)), (Float16(0.2), Float16(0.2)), (Float16(0.3), Float16(0.3))) == Float16(0.0)
            @test orient_fast((0.0, 0.0), (0.1, 0.1), (0.3, 0.2)) < 0.0
            @test orient_fast((0.0, 0.0), (0.1, 0.1), (0.1, 0.2)) > 0.0
            # Case where orient fast fails
            @test orient_fast((0.0, 0.0), (0.1, 0.1), (nextfloat(0.4), 0.4)) > 0.0 # Should be < 0.0
        end
        @testset "Orient Exact" begin
            @test orient_exact((0.1, 0.1), (0.2, 0.2), (0.3, 0.3)) == 0.0
            @test orient_exact((0.0, 3.14159e23), (0.0, 42.4242e11), (0.0, 1.234567e2)) == 0.0
            @test orient_exact((-1.0e-7, -1.0e-7), (2.0, 2.0), (3.0e7, 3.0e7)) == 0.0
            @test orient_exact((Float32(0.1), Float32(0.1)), (Float32(0.2), Float32(0.2)), (Float32(0.3), Float32(0.3))) == Float32(0.0)
            @test orient_exact((Float16(0.1), Float16(0.1)), (Float16(0.2), Float16(0.2)), (Float16(0.3), Float16(0.3))) == Float16(0.0)
            @test orient_exact((0.0, 0.0), (0.1, 0.1), (0.3, 0.2)) < 0.0
            @test orient_exact((0.0, 0.0), (0.1, 0.1), (0.1, 0.2)) > 0.0
            # Case where orient fast fails
            @test orient_exact((0.0, 0.0), (0.1, 0.1), (nextfloat(0.4), 0.4)) < 0.0 # Should be < 0.0
        end
        @testset "Orient" begin
            @test orient((0.1, 0.1), (0.2, 0.2), (0.3, 0.3)) == 0.0
            @test orient((0.0, 3.14159e23), (0.0, 42.4242e11), (0.0, 1.234567e2)) == 0.0
            @test orient((-1.0e-7, -1.0e-7), (2.0, 2.0), (3.0e7, 3.0e7)) == 0.0
            @test orient((Float32(0.1), Float32(0.1)), (Float32(0.2), Float32(0.2)), (Float32(0.3), Float32(0.3))) == Float32(0.0)
            @test orient((Float16(0.1), Float16(0.1)), (Float16(0.2), Float16(0.2)), (Float16(0.3), Float16(0.3))) == Float16(0.0)
            @test orient((0.0, 0.0), (0.1, 0.1), (0.3, 0.2)) < 0.0
            @test orient((0.0, 0.0), (0.1, 0.1), (0.1, 0.2)) > 0.0
            # Case where orient fast fails
            @test orient((0.0, 0.0), (0.1, 0.1), (nextfloat(0.4), 0.4)) < 0.0 # Should be < 0.0
        end
    end

    @testset "Incircle" begin
        @testset "Incircle Fast" begin
            @test incircle_fast((0.5, 0.0), (0.0, 0.5), (-0.5, 0.0), (0.0, 0.0)) > 0.0
            @test incircle_fast((0.5, 0.0), (0.0, 0.5), (-0.5, 0.0), (1.0, 0.0)) < 0.0
            @test incircle_fast((0.5, 0.0), (0.0, 0.5), (-0.5, 0.0), (0.0, -0.5)) == 0.0
            @test incircle_fast((0.5, 0.0), (0.0, 0.5), (-0.5, 0.0), (0.500000001, 0.0)) < 0.0
            @test incircle_fast((0.5, 0.0), (0.0, 0.5), (-0.5, 0.0), (0.499999999, 0.0)) > 0.0
            @test incircle_fast((Float32(0.5), Float32(0.0)), (Float32(0.0), Float32(0.5)), (-Float32(0.5), Float32(0.0)), (Float32(0.0), Float32(0.0))) > Float32(0.0)
            @test incircle_fast((Float32(0.5), Float32(0.0)), (Float32(0.0), Float32(0.5)), (-Float32(0.5), Float32(0.0)), (Float32(1.0), Float32(0.0))) < Float32(0.0)
            @test incircle_fast((Float32(0.5), Float32(0.0)), (Float32(0.0), Float32(0.5)), (-Float32(0.5), Float32(0.0)), (Float32(0.0), -Float32(0.5))) == Float32(0.0)
            @test incircle_fast((Float16(0.5), Float16(0.0)), (Float16(0.0), Float16(0.5)), (-Float16(0.5), Float16(0.0)), (Float16(0.0), Float16(0.0))) > Float16(0.0)
            @test incircle_fast((Float16(0.5), Float16(0.0)), (Float16(0.0), Float16(0.5)), (-Float16(0.5), Float16(0.0)), (Float16(1.0), Float16(0.0))) < Float16(0.0)
            @test incircle_fast((Float16(0.5), Float16(0.0)), (Float16(0.0), Float16(0.5)), (-Float16(0.5), Float16(0.0)), (Float16(0.0), -Float16(0.5))) == Float16(0.0)
            # Degenerate circle
            incircle_fast((0.0, 0.0), (0.5, 0.0), (0.7, 0.0), (0.8, 0.0)) == 0.0
            incircle_fast((0.0, 0.0), (0.5, 0.0), (0.7, 0.0), (0.8, 0.1)) < 0.0
            incircle_fast((0.0, 0.0), (0.5, 0.0), (0.7, 0.0), (0.8, -0.1)) > 0.0
            # Case where incircle fast fails
            @test incircle_fast((0.5, 0.0), (0.0, 0.5), (-0.5, 0.0), (0.0, prevfloat(-0.5))) == 0.0 # should be < 0.0
            @test incircle_fast((0.5, 0.0), (0.0, 0.5), (-0.5, 0.0), (0.0, nextfloat(-0.5))) == 0.0 # should be > 0.0
        end
        @testset "Incircle Exact" begin
            @test incircle_exact((0.5, 0.0), (0.0, 0.5), (-0.5, 0.0), (0.0, 0.0)) > 0.0
            @test incircle_exact((0.5, 0.0), (0.0, 0.5), (-0.5, 0.0), (1.0, 0.0)) < 0.0
            @test incircle_exact((0.5, 0.0), (0.0, 0.5), (-0.5, 0.0), (0.0, -0.5)) == 0.0
            @test incircle_exact((0.5, 0.0), (0.0, 0.5), (-0.5, 0.0), (0.500000001, 0.0)) < 0.0
            @test incircle_exact((0.5, 0.0), (0.0, 0.5), (-0.5, 0.0), (0.499999999, 0.0)) > 0.0
            @test incircle_exact((Float32(0.5), Float32(0.0)), (Float32(0.0), Float32(0.5)), (-Float32(0.5), Float32(0.0)), (Float32(0.0), Float32(0.0))) > Float32(0.0)
            @test incircle_exact((Float32(0.5), Float32(0.0)), (Float32(0.0), Float32(0.5)), (-Float32(0.5), Float32(0.0)), (Float32(1.0), Float32(0.0))) < Float32(0.0)
            @test incircle_exact((Float32(0.5), Float32(0.0)), (Float32(0.0), Float32(0.5)), (-Float32(0.5), Float32(0.0)), (Float32(0.0), -Float32(0.5))) == Float32(0.0)
            @test incircle_exact((Float16(0.5), Float16(0.0)), (Float16(0.0), Float16(0.5)), (-Float16(0.5), Float16(0.0)), (Float16(0.0), Float16(0.0))) > Float16(0.0)
            @test incircle_exact((Float16(0.5), Float16(0.0)), (Float16(0.0), Float16(0.5)), (-Float16(0.5), Float16(0.0)), (Float16(1.0), Float16(0.0))) < Float16(0.0)
            @test incircle_exact((Float16(0.5), Float16(0.0)), (Float16(0.0), Float16(0.5)), (-Float16(0.5), Float16(0.0)), (Float16(0.0), -Float16(0.5))) == Float16(0.0)
            # Degenerate circle
            incircle_exact((0.0, 0.0), (0.5, 0.0), (0.7, 0.0), (0.8, 0.0)) == 0.0
            incircle_exact((0.0, 0.0), (0.5, 0.0), (0.7, 0.0), (0.8, 0.1)) < 0.0
            incircle_exact((0.0, 0.0), (0.5, 0.0), (0.7, 0.0), (0.8, -0.1)) > 0.0
            # Case where incircle fast fails
            @test incircle_exact((0.5, 0.0), (0.0, 0.5), (-0.5, 0.0), (0.0, prevfloat(-0.5))) < 0.0 # should be < 0.0
            @test incircle_exact((0.5, 0.0), (0.0, 0.5), (-0.5, 0.0), (0.0, nextfloat(-0.5))) > 0.0 # should be > 0.0
        end
        @testset "Incircle" begin
            @test incircle((0.5, 0.0), (0.0, 0.5), (-0.5, 0.0), (0.0, 0.0)) > 0.0
            @test incircle((0.5, 0.0), (0.0, 0.5), (-0.5, 0.0), (1.0, 0.0)) < 0.0
            @test incircle((0.5, 0.0), (0.0, 0.5), (-0.5, 0.0), (0.0, -0.5)) == 0.0
            @test incircle((0.5, 0.0), (0.0, 0.5), (-0.5, 0.0), (0.500000001, 0.0)) < 0.0
            @test incircle((0.5, 0.0), (0.0, 0.5), (-0.5, 0.0), (0.499999999, 0.0)) > 0.0
            @test incircle((Float32(0.5), Float32(0.0)), (Float32(0.0), Float32(0.5)), (-Float32(0.5), Float32(0.0)), (Float32(0.0), Float32(0.0))) > Float32(0.0)
            @test incircle((Float32(0.5), Float32(0.0)), (Float32(0.0), Float32(0.5)), (-Float32(0.5), Float32(0.0)), (Float32(1.0), Float32(0.0))) < Float32(0.0)
            @test incircle((Float32(0.5), Float32(0.0)), (Float32(0.0), Float32(0.5)), (-Float32(0.5), Float32(0.0)), (Float32(0.0), -Float32(0.5))) == Float32(0.0)
            @test incircle((Float16(0.5), Float16(0.0)), (Float16(0.0), Float16(0.5)), (-Float16(0.5), Float16(0.0)), (Float16(0.0), Float16(0.0))) > Float16(0.0)
            @test incircle((Float16(0.5), Float16(0.0)), (Float16(0.0), Float16(0.5)), (-Float16(0.5), Float16(0.0)), (Float16(1.0), Float16(0.0))) < Float16(0.0)
            @test incircle((Float16(0.5), Float16(0.0)), (Float16(0.0), Float16(0.5)), (-Float16(0.5), Float16(0.0)), (Float16(0.0), -Float16(0.5))) == Float16(0.0)
            # Degenerate circle
            incircle((0.0, 0.0), (0.5, 0.0), (0.7, 0.0), (0.8, 0.0)) == 0.0
            incircle((0.0, 0.0), (0.5, 0.0), (0.7, 0.0), (0.8, 0.1)) < 0.0
            incircle((0.0, 0.0), (0.5, 0.0), (0.7, 0.0), (0.8, -0.1)) > 0.0
            # Case where incircle fast fails
            @test incircle((0.5, 0.0), (0.0, 0.5), (-0.5, 0.0), (0.0, prevfloat(-0.5))) < 0.0
            @test incircle((0.5, 0.0), (0.0, 0.5), (-0.5, 0.0), (0.0, nextfloat(-0.5))) > 0.0
        end
        @testset "Ghost Incircle" begin
            @testset "Upper Half Plane" begin
                @test incircle((0.42, floatmin()), (0.0, 0.0), (0.5, 0.0)) > 0.0
                @test incircle((0.42, -floatmin()), (0.0, 0.0), (0.5, 0.0)) < 0.0
                @test incircle((floatmin(), 0.42), (0.0, 0.0), (0.0, 0.5)) < 0.0
                @test incircle((-floatmin(), 0.42), (0.0, 0.0), (0.0, 0.5)) > 0.0
            end
            @testset "In Segment" begin
                @test incircle((-floatmin(), 0.0), (0.0, 0.0), (0.5, 0.0)) == 0.0
                @test incircle((0.0, 0.0), (0.0, 0.0), (0.5, 0.0)) == 0.0
                @test incircle((floatmin(), 0.0), (0.0, 0.0), (0.5, 0.0)) == 1.0
                @test incircle((prevfloat(0.5), 0.0), (0.0, 0.0), (0.5, 0.0)) == 1.0
                @test incircle((0.5, 0.0), (0.0, 0.0), (0.5, 0.0)) == 0.0
                @test incircle((nextfloat(0.5), 0.0), (0.0, 0.0), (0.5, 0.0)) == 0.0

                @test incircle((0.0, -floatmin()), (0.0, 0.0), (0.0, 0.5)) == 0.0
                @test incircle((0.0, 0.0), (0.0, 0.0), (0.0, 0.5)) == 0.0
                @test incircle((0.0, floatmin()), (0.0, 0.0), (0.0, 0.5)) == 1.0
                @test incircle((0.0, prevfloat(0.5)), (0.0, 0.0), (0.0, 0.5)) == 1.0
                @test incircle((0.0, 0.5), (0.0, 0.0), (0.0, 0.5)) == 0.0
                @test incircle((0.0, nextfloat(0.5)), (0.0, 0.0), (0.0, 0.5)) == 0.0
            end
        end
    end
end
