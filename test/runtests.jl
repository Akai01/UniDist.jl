using Test
using UniDist

@testset "UniDist.jl" begin
    include("core/test_api.jl")
    include("stats/test_intervals.jl")
    include("continuous/test_continuous.jl")
    include("discrete/test_discrete.jl")
end
