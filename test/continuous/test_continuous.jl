@testset "Continuous distributions" begin
    @testset "StandardUniform" begin
        d = StandardUniform()
        @test pdf(d, 0.5) == 1.0
        @test cdf(d, 0.5) == 0.5
        @test quantile(d, 0.25) == 0.25
        @test support(d) == (0, 1)
    end

    @testset "Exponential" begin
        d = Exponential(2.0)
        @test isapprox(pdf(d, 1.0), 0.5 * exp(-0.5); atol=1e-10)
        @test isapprox(cdf(d, 1.0), 1 - exp(-0.5); atol=1e-10)
        @test isapprox(quantile(d, 0.5), -2.0 * log(0.5); atol=1e-10)
    end

    @testset "Normal" begin
        d = Normal(0.0, 1.0)
        @test isapprox(cdf(d, 0.0), 0.5; atol=1e-10)
        @test isapprox(pdf(d, 0.0), 1 / sqrt(2 * pi); atol=1e-10)
        @test mean(d) == 0.0
        @test var(d) == 1.0
    end

    @testset "Beta" begin
        d = Beta(2.0, 3.0)
        @test isapprox(cdf(d, 0.5), 0.6875; atol=1e-4)
        @test isapprox(quantile(d, 0.5), 0.385; atol=2e-2)
        @test support(d) == (0, 1)
    end

    @testset "Gamma" begin
        d = Gamma(2.0, 2.0)
        @test isapprox(cdf(d, 2.0), 1 - exp(-1) * (1 + 1); atol=2e-3)
    end

    @testset "Uniform" begin
        d = Uniform(0.0, 10.0)
        @test pdf(d, 5.0) == 0.1
        @test cdf(d, 5.0) == 0.5
        @test quantile(d, 0.5) == 5.0
        @test support(d) == (0.0, 10.0)
    end
end
