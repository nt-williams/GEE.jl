module GEE

using LinearAlgebra, Distributions, StatsBase, GLM, Tables
using StatsBase: CoefTable, StatisticalModel, RegressionModel
using LinearAlgebra: BlasReal
using GLM: Link, LinPred

import StatsBase: fit, coef
import GLM: linkfun

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

abstract type GeneralizedEstimatingEquationModel <: StatisticalModel end

include("workingstructures.jl")
include("generalizedestimatingequation.jl")
include("correlation.jl")
include("sandwich.jl")

end