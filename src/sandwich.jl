bread(D, V) = D' * V^-1 * D
meat(D, V, resid) = D' * V^-1 * resid * resid' * V^-1 * D
sandwich(bread::Vector{Matrix{T}}, meat::Vector{Matrix{T}}) where {T<:Real} = sum(bread)^-1 * sum(meat) * sum(bread)^-1

function vcov(m::GeneralizedEstimatingEquationModel)
    V = [updateV(m.R, μ) for μ in m.GEE.μᵢ]
    D = [updateD(m.d, X, μ) for (X, μ) in zip(m.GEE.Xᵢ, m.GEE.μᵢ)]
    rᵢ = [y - μ for (y, μ) in zip(m.GEE.yᵢ, m.GEE.μᵢ)]
    return sandwich(bread.(D, V), meat.(D, V, rᵢ))
end

function stderror(m::GeneralizedEstimatingEquationModel)
    [sqrt(d) for d in diag(vcov(m))]
end