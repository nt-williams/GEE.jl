# # struct GeneralizedEstimatingEquation <: GeneralizedEstimatingEquationsModel
# #     resp
# #     μ̂
# #     β̂
# #     id
# #     corstr
# #     ϕ
# #     vcov
# #     iter
# #     converged
# # end

# # mutable struct GeneralizedEstimatingEquation{G<:GeeResp, L<:LinPred}
# #     rr::G
# #     pp::L
# # end

# """
#     GeeResp
# The response vector and various derived vectors in a GEE model.
# """
# struct GeeResp{V<:AbstractArray{AbstractFloat,1}, D<:UnivariateDistribution, L<:Link}
#     "`y`: response vector"
#     y::V
#     d::D
#     "`eta`: the linear predictor"
#     eta::V
#     "`mu`: mean response"
#     mu::V
# end

# mutable struct GeePred{T<:BlasReal} <: LinPred
#     X::Matrix{T}
#     beta0::Vector{T}
#     iter::Real
#     converged::Bool
# end

# # step 1: generate initial estimates of β̂
# # This needs to fit the GLM and extract the GeneralizedLinearModel object
# function initialfit(::Type{GeneralizedEstimatingEquation}, 
#     f::FormulaTerm, 
#     data,
#     d::UnivariateDistribution = Normal(), 
#     l::Link = canonicallink(d))
#     return glm(f, data, d, l).model
# end

# mutable struct GeeDens
#     yᵢ
#     Xᵢ
#     μᵢ
#     rᵢ
#     β
#     iter
#     converged
# end

# function setup(id, y, X, μ, r, β)
#     uid = unique(id)
#     GeeResp(
#         [y[id .== i] for i in uid], 
#         [X[id .== i, :] for i in uid], 
#         [μ[id .== i] for i in uid], 
#         [r[id .== i] for i in uid], 
#         β, 
#         0, 
#         false
#     )
# end

# # step 2: using current estimates of α (function of working covariance) and ϕ, compute Vᵢ and update β̂
# # slide 34 http://biostat.jhsph.edu/~jleek/teaching/2011/754/lecture7updated.pdf
# function updateβ(β, D, V, r)
#     K = size(D, 1)
#     β + sum([D'[j] * V[j]^-1 * D[j] for j in 1:K])^-1 * sum([D'[j] * V[j]^-1 * r[j] for j in 1:K])
# end

# # step 3: using new β̂, obtain updated α and ϕ using standardized residuals
# # repeat until convergance
# function fit(::Type{GeneralizedEstimatingEquation}, 
#     f::FormulaTerm,
#     data,
#     d::UnivariateDistribution, 
#     l::Link,
#     id, 
#     corstr, 
#     maxit = 20, 
#     tolerance = 0.00001)

#     sch = schema(f, data)
#     f = apply_schema(f, sch)
#     y, X = modelcols(f, data)

#     # step 1: generate initial estimates of β̂
#     initial = initialfit(GeneralizedEstimatingEquation, f, data, d, l)

#     resp = GEE.setup(
#         id,
#         y, 
#         X, 
#         initial.rr.mu, 
#         initial.rr.wrkresid, 
#         initial.pp.beta0
#     )

#     n = length(unique(id))         # total number of cluster
#     m = maximum(size.(resp.Xᵢ, 1)) # number of observations within cluster

#     # step 2: compute α and ϕ
#     # step 3: using current estimates of α (function of working covariance) and ϕ, β̂
#     # Fisher scoring loop, repeat until convergence
#     local ϕ̂, D̂, V̂, R̂
#     while resp.iter <= maxit
#         resp.iter += 1
#         e = pearson.(d, resp.yᵢ, resp.μᵢ)
#         ϕ̂ = ϕ(e)
#         α̂ = updateα(corstr, n, m, e, ϕ̂)
#         R̂ = corstr(m, α̂)
#         V̂ = [V(R̂, μ) for μ in resp.μᵢ]
#         D̂ = [D(d, x, μ) for (x, μ) in zip(resp.Xᵢ, resp.μᵢ)]
#         βnew = GEE.updateβ(resp.β, D̂, V̂, resp.rᵢ)
#         if maximum(abs.((βnew .- resp.β) ./ (resp.β .+ eps(0.0)))) < tolerance
#             resp.converged = true
#             break
#         else
#             resp.β = βnew
#             resp.μᵢ = [X*resp.β for X in resp.Xᵢ]
#             resp.rᵢ = [y - μ for (y, μ) in zip(resp.yᵢ, resp.μᵢ)]
#         end
#     end

#     GeneralizedEstimatingEquation(
#         initial, 
#         collect(Iterators.flatten(resp.μᵢ)), 
#         resp.β, 
#         id, 
#         R̂, 
#         ϕ̂,
#         sandwich(bread.(D̂, V̂), meat.(D̂, V̂, resp.rᵢ)), 
#         resp.iter, 
#         resp.converged
#     )
# end

# coef(m::GeneralizedEstimatingEquation) = m.β

# function StatsBase.coeftable(m::GeneralizedEstimatingEquation)
#     co = coef(m)
#     se = stderror(m)
#     z = co ./ se
#     pvalue = ccdf.(Chisq(1), abs2.(z))
#     CoefTable(
#         hcat(co, se, z, pvalue),
#         ["Coef.", "Std. Error", "z", "Pr(>|z|)"],
#         coefnames(m),
#         4, # pvalcol
#         3, # teststatcol
#     )
# end