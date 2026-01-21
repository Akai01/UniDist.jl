@testset "Statistical intervals" begin
    @testset "Equal-tailed interval" begin
        # Standard normal: 95% interval should be approximately (-1.96, 1.96)
        dist = Normal(0.0, 1.0)
        lo, hi = interval(dist, 0.05)
        @test isapprox(lo, -1.96; atol=0.01)
        @test isapprox(hi, 1.96; atol=0.01)

        # Standard uniform: 90% interval should be (0.05, 0.95)
        dist_u = StandardUniform()
        lo, hi = interval(dist_u, 0.10)
        @test isapprox(lo, 0.05; atol=1e-10)
        @test isapprox(hi, 0.95; atol=1e-10)
    end

    @testset "Highest density interval (continuous)" begin
        # For symmetric distributions, HDI should match equal-tailed interval
        dist = Normal(0.0, 1.0)
        lo_hdi, hi_hdi = hdi(dist, 0.95)
        lo_eq, hi_eq = interval(dist, 0.05)

        # Should be approximately equal for symmetric distributions
        @test isapprox(lo_hdi, lo_eq; atol=0.1)
        @test isapprox(hi_hdi, hi_eq; atol=0.1)
    end

    @testset "Highest density interval (discrete)" begin
        # Binomial distribution
        dist = Binomial(10, 0.5)
        lo, hi = hdi(dist, 0.9)

        # HDI should be within support
        @test lo >= 0
        @test hi <= 10
        @test lo <= hi
    end
end
