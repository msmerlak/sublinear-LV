using DrWatson
@quickactivate


## diversity
function xlogx(x)
    if x <= 0.
        return 0.
    else
        return x*log(x)
    end
end

function diversity(u::Vector)
    v = u./norm(u, 1)
    return exp(-sum(xlogx.(v)))
end
diversity(sol::T) where T<:SciMLBase.AbstractTimeseriesSolution = diversity(sol.u[end])

function richness(u::Vector, tol = 1e-3)
    return length(u[u.>tol])
end
richness(sol::T) where T<:SciMLBase.AbstractTimeseriesSolution = richness(sol.u[end])
