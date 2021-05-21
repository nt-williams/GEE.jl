function StandardizedResiduals(x::GeneralizedEstimatingEquation)
    # (Y_ij - μ_ij) / √(VarianceFunction)
end

VarianceFunction(::Type{Binomial}, μ::Vector{Real}) = μ .* (1 .- μ)

mutable struct V
    μ::Vector{Real}
    corstr::CorrelationStructure
end

(::V) = sqrt(diagm(V.μ)) * V.corstr(length(V.μ))() * sqrt(diagm(μ)) # A^.5 * Corr(Y_ij) * A^0.5

# based on http://biostat.jhsph.edu/~jleek/teaching/2011/754/lecture7updated.pdf
Aᵢ, Bᵢ = [], []
for id in nhefs.seqn
    Xᵢ = init.mm.m[nhefs.seqn .== id, :]
    nᵢ = size(Xᵢ, 1)
    μ̂ᵢ = μ̂[nhefs.seqn .== id]
    rᵢ = r[nhefs.seqn .== id]
    Dᵢ = Xᵢ .* VarianceFunction(Binomial, μ̂ᵢ)
    Vᵢ = sqrt(diagm(μ̂ᵢ)) * IndependenceCorrelation(nᵢ)() * sqrt(diagm(μ̂ᵢ))
    push!(Aᵢ, Dᵢ' * Vᵢ^-1 * Dᵢ)
    push!(Bᵢ, Dᵢ' * Vᵢ^-1 * rᵢ * rᵢ' * Vᵢ^-1 * Dᵢ)
end

sandwich(bread, meat) = bread^-1 * meat * bread^-1
bread(D, V) = D' * V^-1 * D
meat(D, V, resid) = D' * V^-1 * resid * resid' * V^-1 * D

Â = reduce(+, Aᵢ)
B̂ = reduce(+, Bᵢ)

sqrt.(diag(Â^-1 * B̂ * Â^-1))