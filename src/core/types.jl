abstract type AbstractDistribution end
abstract type DiscreteDistribution <: AbstractDistribution end
abstract type ContinuousDistribution <: AbstractDistribution end

support(::AbstractDistribution) = nothing
