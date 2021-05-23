VarianceFunction(::Type{Binomial{T}}, μ̂::Vector{<:Real}) where {T<:Real} = μ̂ .* (1 .- μ̂)

VarianceFunction(::Type{Normal{T}}, μ̂::Vector{<:Real}) where {T<:Real} = 1

VarianceFunction(::Type{Poisson{T}}, μ̂::Vector{<:Real}) where {T<:Real} = μ̂

pearson(d::UnivariateDistribution, Y, μ̂) = (Y .- μ̂) ./ sqrt(VarianceFunction(typeof(d), μ̂))

function ϕ(e)
    phi, m, n = 0, length(e), maximum(length.(e))
    for i in 1:m
        for j in 1:length(e[i])
            phi += e[i][j]^2
        end
    end
    phi / (n * m)
end

# D' = X'A where A is diagonal matrix with Var(Y_ij|X_ij), i.e., variance function
D(d::UnivariateDistribution, X::Matrix, μ̂::Vector) = X .* VarianceFunction(typeof(d), μ̂)

V(corstr::CorrelationStructure, μ̂::Vector) = sqrt(diagm(μ̂)) * corstr() * sqrt(diagm(μ̂))