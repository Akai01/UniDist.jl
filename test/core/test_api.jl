@testset "Core API" begin
    @testset "Shorthand functions" begin
        dist = StandardUniform()

        # d/p/q shortcuts should match pdf/cdf/quantile
        @test d(dist, 0.5) == pdf(dist, 0.5)
        @test p(dist, 0.5) == cdf(dist, 0.5)
        @test q(dist, 0.5) == quantile(dist, 0.5)
    end

    @testset "Survival and hazard functions" begin
        dist = Exponential(1.0)
        x = 1.0

        # Survival function: sf(x) = 1 - cdf(x)
        @test isapprox(sf(dist, x), 1 - cdf(dist, x); atol=1e-10)

        # Hazard function: hazard(x) = pdf(x) / sf(x)
        @test isapprox(hazard(dist, x), pdf(dist, x) / sf(dist, x); atol=1e-10)

        # Cumulative hazard: cumhaz(x) = -log(sf(x))
        @test isapprox(cumhaz(dist, x), -log(sf(dist, x)); atol=1e-10)
    end

    @testset "Random sampling" begin
        dist = Normal(0.0, 1.0)

        # r() should return array of samples
        samples = r(dist, 100)
        @test length(samples) == 100
        @test eltype(samples) <: Real
    end
end
