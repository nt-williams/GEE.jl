abstract type CorrelationStructure end

# not mutable because it will never change
struct Independence <: CorrelationStructure
    N::Number
    α
end

Independence(N, α = 0) = Independence(N, α)

mutable struct Exchangeable <: CorrelationStructure
    N::Real
    α̂::Real
end

mutable struct Unstructered <: CorrelationStructure
    N::Real
    α̂
end

(X::Independence)() = Matrix{Int}(I, X.N, X.N)

(X::Unstructered)() = reshape(X.α̂, (X.N, X.N))

function (X::Exchangeable)()
    R = zeros(X.N, X.N)
    for j in 1:X.N
        for k in 1:X.N
            R[j, k] = j == k ? 1 : X.α̂
        end
    end
    return R
end

function updateα(::Type{Independence}, n, m, e, ϕ)
    Independence(m)()
end

function updateα(::Type{Exchangeable}, n, m, e, ϕ)
    R = updateα(Unstructered, n, m, e, ϕ)
    mean([R[r] for r in CartesianIndices(R) if r[1] ≠ r[2]])
end

# n is the number of clusters
# m is the number of observations within clusters
# e is the pearsons residuals within a cluster
# ϕ is well... ϕ
# returns the unstructured correlation matrix
function updateα(::Type{Unstructered}, n, m, e, ϕ)
    R = zeros(m, m)
    for i in 1:n
        for j in 1:m
            for k in 1:m
                R[j, k] += e[i][j]*e[i][k]
            end
        end
    end
    R = R / (n * ϕ)
    R[diagind(R)] .= 1
    return R
end
