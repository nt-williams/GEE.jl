bread(D, V) = D' * V^-1 * D

meat(D, V, resid) = D' * V^-1 * resid * resid' * V^-1 * D

sandwich(bread::Vector{Matrix{T}}, meat::Vector{Matrix{T}}) where {T<:Real} = sum(bread)^-1 * sum(meat) * sum(bread)^-1

robustse(vcov) = sqrt.(diag(vcov))
