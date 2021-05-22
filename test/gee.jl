include("src/GEE.jl")

using .GEE
using GLM
using RCall

R"""
library(gee)
data(warpbreaks)

fit <- summary(gee(breaks ~ tension, id=wool, data=warpbreaks, corstr="exchangeable"))

library(geeM)
summary(geem(breaks ~ tension, id=wool, data=warpbreaks, corstr="exchangeable", useP = FALSE))
"""

warpbreaks = rcopy(R"warpbreaks")

fitted = GEE.initialfit(
    GEE.GeneralizedEstimatingEquation, 
    @formula(breaks ~ tension), 
    warpbreaks,
    Normal(), 
    IdentityLink()
)

id = warpbreaks.wool
Yᵢ = [warpbreaks.breaks[id .== i, :] for i in unique(id)]
Xᵢ = [fitted.pp.X[id .== i, :] for i in unique(id)]
μ̂ᵢ = [predict(fitted)[id .== i] for i in unique(id)]
residᵢ = [warpbreaks.breaks[id .== i] .- predict(fitted)[id .== i] for i in unique(id)]

# this should be called Pearsons
e = GEE.StandardizedResiduals.(Normal(), Yᵢ, μ̂ᵢ)

N = length(id)
n = length(unique(id))
m = length(Yᵢ[1])
r = zeros(m, m)

phi = GEE.ϕ(e)
# phi = 141.1481 # when useP = TRUE in geeM

# this will create the unstructured correlation matrix
function alpha(n, m, e, phi)
    R = zeros(m, m)
    for i in 1:n
        for j in 1:m
            for k in 1:m
                if j == k
                    R[j, k] = 1
                else
                    R[j, k] += e[i][j]*e[i][k]
                end
            end
        end
    end
    R = R / (n * phi)
    R[diagind(R)] .= 1
    return R
end

unstr = alpha(n, m, e, phi)

# for i in 1:n
#     for j in 1:27
#         for k in 1:27
#             if j == k
#                 r[j, k] = 1
#             else
#                 r[j, k] += e[i][j]*e[i][k]
#             end
#         end
#     end
# end

# unstr = r / (n * phi)

function offdiag2(A::AbstractMatrix)
    [A[ι] for ι in CartesianIndices(A) if ι[1] ≠ ι[2]]
end 

# we can calculate the exchangeable alpha from the unstructured correlation
a = mean(offdiag2(r / 2phi)) # I am successfully calculating exchangeable alpha when useP = FALSE in geeM

# R should theoretically be the correlation matrix, R()
# R = GEE.Independence(m)
R = GEE.Exchangeable(m, a)

Vᵢ = [GEE.V(R, mu) for mu in μ̂ᵢ]
Dᵢ = GEE.D.(Normal(), Xᵢ, μ̂ᵢ)

gee = GEE.GeeComponents(n, fitted.pp.beta0, Dᵢ, Vᵢ, residᵢ, 0, 0, 0)
GEE.updateβ(gee)

A = GEE.bread.(Dᵢ, Vᵢ)
B = GEE.meat.(Dᵢ, Vᵢ, residᵢ)
GEE.sandwich(A, B)
GEE.robustse(A, B)