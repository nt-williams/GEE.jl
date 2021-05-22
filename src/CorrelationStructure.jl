abstract type CorrelationStructure end

# shouldn't be mutable because it will never change
struct IndependenceCorrelation <: CorrelationStructure
    N::Number
end

mutable struct Exchangeable <: CorrelationStructure
    N::Real
    α̂::Real
end

mutable struct Unstructered <: CorrelationStructure
    N::Real
end

(X::IndependenceCorrelation)() = Matrix{Int}(I, X.N, X.N)

function (X::Exchangeable)()
    R = zeros(X.N, X.N)
    for j in 1:X.N
        for k in 1:X.N
            if j == k
                R[j, k] = 1
            else
                R[j, k] = X.α̂
            end
        end
    end
    return R
end

function updateα(::Type{Exchangeable}, e, ϕ)
end

# function updateα(::Type{Unstructered}, gee::GeeComponents)
# end
