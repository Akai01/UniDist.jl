@testset "Discrete distributions" begin
    @testset "Bernoulli" begin
        d = Bernoulli(0.3)
        @test pdf(d, 0) == 0.7
        @test pdf(d, 1) == 0.3
        @test cdf(d, 0) == 0.7
        @test cdf(d, 1) == 1.0
        @test mean(d) == 0.3
    end

    @testset "Binomial" begin
        d = Binomial(10, 0.5)
        @test isapprox(pdf(d, 5), 0.24609375; atol=1e-6)
        @test isapprox(mean(d), 5.0; atol=1e-10)
        @test isapprox(var(d), 2.5; atol=1e-10)
    end

    @testset "Poisson" begin
        d = Poisson(3.0)
        # P(X=0) = e^(-3)
        @test isapprox(pdf(d, 0), exp(-3); atol=1e-10)
        # P(X=3) = e^(-3) * 3^3 / 3!
        @test isapprox(pdf(d, 3), exp(-3) * 27 / 6; atol=1e-10)
        @test mean(d) == 3.0
        @test var(d) == 3.0
    end

    @testset "DiscreteUniform" begin
        d = DiscreteUniform(1, 6)
        @test isapprox(pdf(d, 3), 1/6; atol=1e-10)
        @test isapprox(mean(d), 3.5; atol=1e-10)
    end

    @testset "Geometric" begin
        d = Geometric(0.5)
        # P(X=0) = 0.5 (first trial success)
        @test pdf(d, 0) == 0.5
        # P(X=1) = 0.5 * 0.5 = 0.25
        @test pdf(d, 1) == 0.25
    end
end
