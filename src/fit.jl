abstract type GeneralizedEstimatingEquationsModel{T} <: StatsModels.RegressionModel end

mutable struct GeneralizedEstimatingEquation{T} <: GeneralizedEstimatingEquationsModel{T}
    resp::GLM.GlmResp
    corstr::CorrelationStructure
end

# step 1: generate initial estimates of β̂
function initialfit(::Type{GeneralizedEstimatingEquation}, 
    f::FormulaTerm, 
    data,
    d::UnivariateDistribution = Normal(), 
    l::Link = canonicallink(d))
    # This needs to fit the GLM and extract the GlmResp object
end

# step 2: using current estimates of α (function of working covariance) and ϕ, compute Vᵢ and update β̂
function updateβ!(x::GeneralizedEstimatingEquation)
    # solve for a new β by using equation 13.2 (pg 356) in https://onlinelibrary.wiley.com/doi/pdf/10.1002/9781119513469.ch13
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