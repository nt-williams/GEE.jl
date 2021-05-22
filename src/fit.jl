abstract type GeneralizedEstimatingEquationsModel{T} <: StatsModels.RegressionModel end

mutable struct GeneralizedEstimatingEquation{T} <: GeneralizedEstimatingEquationsModel{T}
    resp::GLM.GlmResp
    corstr::CorrelationStructure
    id::Vector
    y::Vector{<:Real}
    μ̂::Vector{<:Real}
    β̂::Vector{<:Real}
end

# step 1: generate initial estimates of β̂
function initialfit(::Type{GeneralizedEstimatingEquation}, 
    f::FormulaTerm, 
    data,
    d::UnivariateDistribution = Normal(), 
    l::Link = canonicallink(d))
    return glm(f, data, d, l).model
    # This needs to fit the GLM and extract the GeneralizedLinearModel object
end

mutable struct GeeComponents
    K
    β̂
    D
    V
    resid
    e
    ϕ
    iter
    converged
end

# step 2: using current estimates of α (function of working covariance) and ϕ, compute Vᵢ and update β̂
# slide 34 http://biostat.jhsph.edu/~jleek/teaching/2011/754/lecture7updated.pdf
function updateβ(gee::GeeComponents)
    gee.β̂ + sum([gee.D'[j] * gee.V[j]^-1 * gee.D[j] for j in 1:gee.K])^-1 * sum([gee.D'[j] * gee.V[j]^-1 * gee.resid[j] for j in 1:gee.K])
end

# step 3: using new β̂, obtain updated α and ϕ using standardized residuals
# repeat until convergance

# function fit(::Type{GeneralizedEstimatingEquation}, 
#     f::FormulaTerm,
#     data,
#     d::UnivariateDistribution = Normal(), 
#     l::Link = canonicallink(d),
#     id::AbstractVector{<:Real}, 
#     corstr::CorrelationStructure = IndependenceCorrelation;
#     dofit::Bool = true, 
#     wts::AbstractVector{<:Real} = similar(y, 0),
#     offset::AbstractVector{<:Real} = similar(y, 0))
# end