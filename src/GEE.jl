module GEE

using LinearAlgebra
using GLM
using StatsModels
using Distributions

using GLM: Link

export GeneralizedEstimatingEquationsModel
export GeneralizedEstimatingEquation
export CorrelationStructure
export Independence
export Exchangeable
export Unstructered
export fit

include("workingstructures.jl")
include("fit.jl")
include("correlation.jl")
include("sandwich.jl")

end