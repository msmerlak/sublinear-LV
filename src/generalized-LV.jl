using DrWatson
@quickactivate

include(srcdir("utils.jl"))

using DifferentialEquations
using Random, Distributions
using LinearAlgebra

using SparseArrays

## generalized LV dynamics

function interaction_matrix(S, μ, σ, scaled = true, symm = false)
    scaled ? A = rand(Normal(μ/S, σ/sqrt(S)), (S,S)) : A = rand(Normal(μ, σ), (S,S))
    if symm
        A = (A + A')/sqrt(2)
    end
    A[diagind(A)] .= 0.
    return A
end

function F(u, p)
    return p["r"].*(u.^p["k"] .- u.^2. / p["K"]^(2-p["k"])) .- u.*(p["A"]*u) 
end


function drift!(du, u, p, t)
        u .= ppart.(u)
        du .= p["r"].*(u.^p["k"] .- u.^2. / p["K"]^(2-p["k"])) .- u.*(p["A"]*u) .+ p["imm"]
end


function diffusion!(du, u, p, t)
    u .= ppart.(u)
    du .= sqrt.(p["T"]*u)
end


## solving
blowup(threshold = 500) = DiscreteCallback((u, t, integrator) -> maximum(u) > threshold, terminate!)
collapse(fraction = 0.1) = DiscreteCallback((u, t, integrator) -> diversity(u) < fraction*length(u), terminate!)
converged(abstol = 1e-6) = TerminateSteadyState(abstol)


function  LV_problem(params; k = 1., x₀ = nothing, max_time = 1000.)
    
    p = copy(params)
    p["imm"] = fill(p["λ"], p["S"])
    if "A" ∉ keys(p)
        p["A"] = interaction_matrix(p["S"], p["μ"], p["σ"], p["scaled"], p["symm"])
    end

    @assert p["ξ"] < 1 "Variance of growth rates too large"
    p["r"] = 1 .+ p["ξ"]*(rand(p["S"]) .- 0.5)

    if x₀ === nothing
        x₀ = rand(p["S"]).* (p["K"]/2)
    end
    
    p["k"] = k

    if p["T"] == 0.
        pb = ODEProblem(drift!, x₀, (0., max_time), p);
    else
        pb = SDEProblem(drift!, diffusion!, x₀, (0., max_time), p);
    end
    
    return (pb, p)
end

function  LV_solver(params; k = 1., x₀ = nothing, termination = [blowup()], max_time = 1000.)

    pb, params = LV_problem(params, k = k, x₀ = x₀, max_time = max_time)
    return solve(pb,
    saveat = 0:1:max_time,
     callback = CallbackSet(termination...)), params
end

## lyapunov spectrum

using DynamicalSystems
function  LV_lyap(params; k = 1, x₀ = nothing, termination = [blowup()], max_time = 1000.)


    @assert p["T"] == 0. "Lyapunov spectrum only at zero T"
    
    pb, params = LV_problem(params, k = k, x₀ = x₀, max_time = max_time)

    return lyapunov(ContinuousDynamicalSystem(pb), max_time, 
    alg = BS3(),
    callback = CallbackSet(termination...))
end