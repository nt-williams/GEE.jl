module GEE

using LinearAlgebra
using GLM
using StatsModels
using StatsBase
using Distributions

using GLM: Link

import StatsBase: fit, coef

export GeneralizedEstimatingEquationsModel,
       GeneralizedEstimatingEquation,
       CorrelationStructure,
       Independence,
       Exchangeable,
       Unstructered,
       fit,
       vcov,
       coef,
       stderror,
       Normal, 
       Binomial, 
       Poisson,
       LogitLink, 
       IdentityLink, 
       LogLink,
       @formula

include("workingstructures.jl")
include("fit.jl")
include("correlation.jl")
include("sandwich.jl")

end