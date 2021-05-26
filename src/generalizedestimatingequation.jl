struct GeneralizedEstimatingEquation <: GeneralizedEstimatingEquationModel
    GEE
    β
    ϕ::Real
    μ
    η
    d::Distribution
    l::Link
    R
end

struct Model{T<:AbstractFloat}
    X::Matrix{T}
    y::Vector{T}
    d::Distribution
    l::Link
end

mutable struct GeeDens
    yᵢ::Vector
    Xᵢ::Vector
    μᵢ::Vector
    rᵢ::Vector
    iter::Real
    converged::Bool
end

function GeeModel(f::FormulaTerm, data, d::Distribution, l::Link, contrasts = Dict{Symbol,Any}())
    sch = schema(f, data, contrasts)
    f = apply_schema(f, sch)
    y, X = modelcols(f, data)
    Model(X, y, d, l)
end

function GeeDens(y, X, μ, r, id)
    uid = unique(id)
    Xᵢ, yᵢ, μᵢ, rᵢ = [], [], [], []
    for ident in uid
        keep = id .== ident
        push!(Xᵢ, X[keep, :])
        push!(yᵢ, y[keep])
        push!(μᵢ, μ[keep])
        push!(rᵢ, r[keep])
    end
    GeeDens(yᵢ, Xᵢ, μᵢ, rᵢ, 0, false)
end

Base.size(x::GeeDens) = (n = size(x.Xᵢ, 1), m =  maximum(size.(x.Xᵢ, 1)))

function GeneralizedEstimatingEquation(    
    f::FormulaTerm,
    data,
    d::Distribution,
    l::Link = canonicallink(d), 
    id = [], 
    corstr::CorrelationStructure = Independence;
    contrasts = Dict{Symbol,Any}(), 
    maxitr::Integer = 10, 
    tol::Real = 0.0001
    )

    Xy = GeeModel(f, data, d, l, contrasts)  # extract X and y from a formula
    ifit = glm(Xy.X, Xy.y, d, l)            # conduct initial fit using glm
    β = coef(ifit)                          # extract initial coefficient estimates
    params = GeeDens(Xy.X, Xy.y, ifit.rr.mu, ifit.rr.wrkresid, id) # create Dens predict object where necessary components split by id

    GeneralizedEstimatingEquation(params, d, l, β, corstr, maxitr = maxitr, tol = tol)
end

function GeneralizedEstimatingEquation(
    x, 
    d::Distribution, 
    l::Link,
    beta, 
    corstr; 
    maxitr::Integer = 10, 
    tol::Real = 0.0001)

    β = copy(beta)
    n, m = size(x)

    local ϕ, D, V, R
    while x.iter <= maxitr
        x.iter += 1

        e = pearson.(d, x.yᵢ, x.μᵢ)
        ϕ = phi(e)
        α = updateα(corstr, n, m, e, ϕ)
        R = corstr(m, α)
        V = [updateV(R, μ) for μ in x.μᵢ]
        D = [updateD(d, X, μ) for (X, μ) in zip(x.Xᵢ, x.μᵢ)]
        βnew = updateβ(β, D, V, x.rᵢ)

        if maximum(abs.((βnew .- β) ./ (β .+ eps(0.0)))) < tol
            x.converged = true
            break
        else
            β = βnew
            x.μᵢ = [X*β for X in x.Xᵢ]
            x.rᵢ = [y - μ for (y, μ) in zip(x.yᵢ, x.μᵢ)]
        end
    end

    μ = collect(Iterators.flatten(x.μᵢ))
    η = linkfun.(l, μ)
    return GeneralizedEstimatingEquation(x, β, ϕ, μ, η, d, l, R)
end

function updateβ(β, D, V, r)
    K = size(D, 1)
    β + sum([D'[j] * V[j]^-1 * D[j] for j in 1:K])^-1 * sum([D'[j] * V[j]^-1 * r[j] for j in 1:K])
end