abstract type CorrelationStructure end

struct IndependenceCorrelation <: CorrelationStructure
    N::Number
end

(X::IndependenceCorrelation)() = diagm([1 for n in 1:X.N])
