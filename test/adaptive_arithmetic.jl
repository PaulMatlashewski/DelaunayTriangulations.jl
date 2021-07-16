using DelaunayTriangulations.GeometryPredicates:
    two_sum_tail,
    two_sum,
    fast_two_sum_tail,
    fast_two_sum,
    two_diff_tail,
    two_diff,
    two_prod_tail,
    two_prod,
    grow_expansion,
    expansion_sum,
    expansion_diff,
    expansion_scale,
    approximate

@testset "Adaptive Arithmetic" begin
    @testset "Two Sum" begin
        @testset "Float64" begin
            @test two_sum_tail(0.1, 0.2, 0.1 + 0.2) < 0.0
            @test two_sum_tail(0.2, 0.2, 0.2 + 0.2) == 0.0
            @test two_sum_tail(0.4, 0.3, 0.4 + 0.3) > 0.0
            @test two_sum_tail(-0.1, -0.2, -0.1 - 0.2) > 0.0
            @test two_sum(0.1, 0.2) == (two_sum_tail(0.1, 0.2, 0.1 + 0.2), 0.1 + 0.2)
            @test two_sum(0.2, 0.2) == (two_sum_tail(0.2, 0.2, 0.2 + 0.2), 0.2 + 0.2)
            @test two_sum(0.4, 0.3) == (two_sum_tail(0.4, 0.3, 0.4 + 0.3), 0.4 + 0.3)
            @test two_sum(-0.1, -0.2) == (two_sum_tail(-0.1, -0.2, -0.1 - 0.2), -0.1 - 0.2)
        end
        @testset "Float32" begin
            @test two_sum_tail(Float32(0.1), Float32(0.2), Float32(0.1) + Float32(0.2)) < Float32(0.0)
            @test two_sum_tail(Float32(0.2), Float32(0.2), Float32(0.2) + Float32(0.2)) == Float32(0.0)
            @test two_sum_tail(Float32(0.3), Float32(0.2), Float32(0.3) + Float32(0.2)) > Float32(0.0)
            @test two_sum_tail(-Float32(0.1), -Float32(0.2), -Float32(0.1) - Float32(0.2)) > Float32(0.0)
            @test two_sum(Float32(0.1), Float32(0.2)) == (two_sum_tail(Float32(0.1), Float32(0.2), Float32(0.1) + Float32(0.2)), Float32(0.1) + Float32(0.2))
            @test two_sum(Float32(0.2), Float32(0.2)) == (two_sum_tail(Float32(0.2), Float32(0.2), Float32(0.2) + Float32(0.2)), Float32(0.2) + Float32(0.2))
            @test two_sum(Float32(0.3), Float32(0.2)) == (two_sum_tail(Float32(0.3), Float32(0.2), Float32(0.3) + Float32(0.2)), Float32(0.3) + Float32(0.2))
            @test two_sum(-Float32(0.1), -Float32(0.2)) == (two_sum_tail(-Float32(0.1), -Float32(0.2), -Float32(0.1) - Float32(0.2)), -Float32(0.1) - Float32(0.2))
        end
        @testset "Float16" begin
            @test two_sum_tail(Float16(0.1), Float16(0.2), Float16(0.1) + Float16(0.2)) > Float16(0.0)
            @test two_sum_tail(Float16(0.2), Float16(0.2), Float16(0.2) + Float16(0.2)) == Float16(0.0)
            @test two_sum_tail(Float16(0.5), Float16(0.2), Float16(0.5) + Float16(0.2)) < Float16(0.0)
            @test two_sum_tail(-Float16(0.1), -Float16(0.2), -Float16(0.1) - Float16(0.2)) < Float16(0.0)
            @test two_sum(Float16(0.1), Float16(0.2)) == (two_sum_tail(Float16(0.1), Float16(0.2), Float16(0.1) + Float16(0.2)), Float16(0.1) + Float16(0.2))
            @test two_sum(Float16(0.2), Float16(0.2)) == (two_sum_tail(Float16(0.2), Float16(0.2), Float16(0.2) + Float16(0.2)), Float16(0.2) + Float16(0.2))
            @test two_sum(Float16(0.5), Float16(0.2)) == (two_sum_tail(Float16(0.5), Float16(0.2), Float16(0.5) + Float16(0.2)), Float16(0.5) + Float16(0.2))
            @test two_sum(-Float16(0.1), -Float16(0.2)) == (two_sum_tail(-Float16(0.1), -Float16(0.2), -Float16(0.1) - Float16(0.2)), -Float16(0.1) - Float16(0.2))
        end
    end
    @testset "Fast Two Sum" begin
        @testset "Float64" begin
            @test fast_two_sum_tail(0.2, 0.1, 0.2 + 0.1) < 0.0
            @test fast_two_sum_tail(0.4, 0.3, 0.4 + 0.3) > 0.0
            @test fast_two_sum_tail(-0.1, -0.2, -0.1 - 0.2) > 0.0
            @test fast_two_sum(0.2, 0.1) == (fast_two_sum_tail(0.2, 0.1, 0.2 + 0.1), 0.2 + 0.1)
            @test fast_two_sum(0.4, 0.3) == (fast_two_sum_tail(0.4, 0.3, 0.4 + 0.3), 0.4 + 0.3)
            @test fast_two_sum(-0.1, -0.2) == (fast_two_sum_tail(-0.1, -0.2, -0.1 - 0.2), -0.1 - 0.2)
        end
        @testset "Float32" begin
            @test fast_two_sum_tail(Float32(0.2), Float32(0.1), Float32(0.2) + Float32(0.1)) < Float32(0.0)
            @test fast_two_sum_tail(Float32(0.3), Float32(0.2), Float32(0.3) + Float32(0.2)) > Float32(0.0)
            @test fast_two_sum_tail(-Float32(0.1), -Float32(0.2), -Float32(0.1) - Float32(0.2)) > Float32(0.0)
            @test fast_two_sum(Float32(0.2), Float32(0.1)) == (fast_two_sum_tail(Float32(0.2), Float32(0.1), Float32(0.2) + Float32(0.1)), Float32(0.2) + Float32(0.1))
            @test fast_two_sum(Float32(0.3), Float32(0.2)) == (fast_two_sum_tail(Float32(0.3), Float32(0.2), Float32(0.3) + Float32(0.2)), Float32(0.3) + Float32(0.2))
            @test fast_two_sum(-Float32(0.1), -Float32(0.2)) == (fast_two_sum_tail(-Float32(0.1), -Float32(0.2), -Float32(0.1) - Float32(0.2)), -Float32(0.1) - Float32(0.2))
        end
        @testset "Float16" begin
            @test fast_two_sum_tail(Float16(0.2), Float16(0.1), Float16(0.2) + Float16(0.1)) > Float16(0.0)
            @test fast_two_sum_tail(Float16(0.5), Float16(0.2), Float16(0.5) + Float16(0.2)) < Float16(0.0)
            @test fast_two_sum_tail(-Float16(0.1), -Float16(0.2), -Float16(0.1) - Float16(0.2)) < Float16(0.0)
            @test fast_two_sum(Float16(0.2), Float16(0.1)) == (fast_two_sum_tail(Float16(0.2), Float16(0.1), Float16(0.2) + Float16(0.1)), Float16(0.2) + Float16(0.1))
            @test fast_two_sum(Float16(0.5), Float16(0.2)) == (fast_two_sum_tail(Float16(0.5), Float16(0.2), Float16(0.5) + Float16(0.2)), Float16(0.5) + Float16(0.2))
            @test fast_two_sum(-Float16(0.1), -Float16(0.2)) == (fast_two_sum_tail(-Float16(0.1), -Float16(0.2), -Float16(0.1) - Float16(0.2)), -Float16(0.1) - Float16(0.2))
        end
    end
    @testset "Two Diff" begin
        @testset "Float64" begin
            @test two_diff_tail(0.3, 0.11, 0.3 - 0.11) < 0.0
            @test two_diff_tail(0.3, 0.1, 0.3 - 0.1) == 0.0
            @test two_diff_tail(0.411, 0.92, 0.411 - 0.92) > 0.0
            @test two_diff(0.3, 0.11) == (two_diff_tail(0.3, 0.11, 0.3 - 0.11), 0.3 - 0.11)
            @test two_diff(0.3, 0.1) == (two_diff_tail(0.3, 0.1, 0.3 - 0.1), 0.3 - 0.1)
            @test two_diff(0.411, 0.92) == (two_diff_tail(0.411, 0.92, 0.411 - 0.92), 0.411 - 0.92)
        end
        @testset "Float32" begin
            @test two_diff_tail(Float32(0.3), Float32(0.1), Float32(0.3) - Float32(0.1)) < Float32(0.0)
            @test two_diff_tail(Float32(0.3), Float32(0.11), Float32(0.3) - Float32(0.11)) == Float32(0.0)
            @test two_diff(Float32(0.3), Float32(0.1)) == (two_diff_tail(Float32(0.3), Float32(0.1), Float32(0.3) - Float32(0.1)), Float32(0.3) - Float32(0.1))
            @test two_diff(Float32(0.3), Float32(0.11)) == (two_diff_tail(Float32(0.3), Float32(0.11), Float32(0.3) - Float32(0.11)), Float32(0.3) - Float32(0.11))
        end
        @testset "Float16" begin
            @test two_diff_tail(Float16(0.11), Float16(0.412), Float16(0.11) - Float16(0.412)) > Float16(0.0)
            @test two_diff_tail(Float16(0.3), Float16(0.11), Float16(0.3) - Float16(0.11)) == Float16(0.0)
            @test two_diff(Float16(0.3), Float16(0.1)) == (two_diff_tail(Float16(0.3), Float16(0.1), Float16(0.3) - Float16(0.1)), Float16(0.3) - Float16(0.1))
            @test two_diff(Float16(0.3), Float16(0.11)) == (two_diff_tail(Float16(0.3), Float16(0.11), Float16(0.3) - Float16(0.11)), Float16(0.3) - Float16(0.11))
        end
    end
    @testset "Two Prod" begin
        @testset "Float64" begin
            @test two_prod_tail(0.3, 0.11, 0.3 * 0.11) < 0.0
            @test two_prod_tail(0.3, 0.1, 0.3 * 0.1) > 0.0
            @test two_prod_tail(2.0, 2.0, 2.0 * 2.0) == 0.0
            @test two_prod(0.3, 0.11) == (two_prod_tail(0.3, 0.11, 0.3 * 0.11), 0.3 * 0.11)
            @test two_prod(0.3, 0.1) == (two_prod_tail(0.3, 0.1, 0.3 * 0.1), 0.3 * 0.1)
            @test two_prod(0.2, 0.2) == (two_prod_tail(0.2, 0.2, 0.2 * 0.2), 0.2 * 0.2)
        end
        @testset "Float32" begin
            @test two_prod_tail(Float32(0.3), Float32(0.11), Float32(0.3) * Float32(0.11)) > Float32(0.0)
            @test two_prod_tail(Float32(0.3), Float32(0.09), Float32(0.3) * Float32(0.09)) < Float32(0.0)
            @test two_prod_tail(Float32(2.0), Float32(2.0), Float32(2.0) * Float32(2.0)) == Float32(0.0)
            @test two_prod(Float32(0.3), Float32(0.11)) == (two_prod_tail(Float32(0.3), Float32(0.11), Float32(0.3) * Float32(0.11)), Float32(0.3) * Float32(0.11))
            @test two_prod(Float32(0.3), Float32(0.09)) == (two_prod_tail(Float32(0.3), Float32(0.09), Float32(0.3) * Float32(0.09)), Float32(0.3) * Float32(0.09))
            @test two_prod(Float32(0.2), Float32(0.2)) == (two_prod_tail(Float32(0.2), Float32(0.2), Float32(0.2) * Float32(0.2)), Float32(0.2) * Float32(0.2))
        end
        @testset "Float16" begin
            @test two_prod_tail(Float16(0.3), Float16(0.11), Float16(0.3) * Float16(0.11)) > Float16(0.0)
            @test two_prod_tail(Float16(0.3), Float16(0.05), Float16(0.3) * Float16(0.05)) < Float16(0.0)
            @test two_prod_tail(Float16(2.0), Float16(2.0), Float16(2.0) * Float16(2.0)) == Float16(0.0)
            @test two_prod(Float16(0.3), Float16(0.11)) == (two_prod_tail(Float16(0.3), Float16(0.11), Float16(0.3) * Float16(0.11)), Float16(0.3) * Float16(0.11))
            @test two_prod(Float16(0.3), Float16(0.05)) == (two_prod_tail(Float16(0.3), Float16(0.05), Float16(0.3) * Float16(0.05)), Float16(0.3) * Float16(0.05))
            @test two_prod(Float16(0.2), Float16(0.2)) == (two_prod_tail(Float16(0.2), Float16(0.2), Float16(0.2) * Float16(0.2)), Float16(0.2) * Float16(0.2))
        end
    end
    @testset "Expansions" begin
        @testset "Grow Expansion" begin
            @test grow_expansion(0.123, 0.456) == two_sum(0.123, 0.456)
            @test grow_expansion(0.42) == 0.42
            e = rand()
            for i in 1:10
                e = grow_expansion(e, rand())
                @test length(e) == i + 1
                for j in 1:i
                    @test e[j] < e[j+1] || e[j] == 0.0 || e[j+1] == 0.0
                end
            end
        end
        @testset "Expansion Sum" begin
            for i in 1:3
                e = grow_expansion(([rand() for _ in 1:i]...,), rand())
                for j in 1:3
                    f = grow_expansion(([rand() for _ in 1:j]...,), rand())
                    g = expansion_sum(e, f)
                    @test length(g) == length(e) + length(f)
                    @test g[end] ≈ e[end] + f[end]
                end
            end
        end
        @testset "Expansion Diff" begin
            e, f = grow_expansion(rand(), rand()), rand()
            g = expansion_diff(e, f)
            @test length(g) == length(e) + length(f)
            @test g[end] ≈ e[end] - f[end]
            e, f = grow_expansion(rand(), rand()), grow_expansion(rand(), rand())
            g = expansion_diff(e, f)
            @test length(g) == length(e) + length(f)
            @test g[end] ≈ e[end] - f[end]
        end
        @testset "Expansion Scale" begin
            for i in 1:10
                e, f = grow_expansion(([rand() for _ in 1:i]...,), rand()), rand()
                g = expansion_scale(e, f)
                @test length(g) == 2 * length(e)
                @test sum(g) ≈ f * sum(e)
            end
        end
        @testset "Approximate" begin
            e = grow_expansion(([rand() for _ in 1:10]...,), rand())
            @test approximate(e) == sum(e)
        end
    end
end
