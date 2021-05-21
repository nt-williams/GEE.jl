module GEE

using LinearAlgebra
using GLM
using StatsModels
using Distributions

using GLM: Link

export IndependenceCorrelation

include("CorrelationStructure.jl")
include("fit.jl")

end