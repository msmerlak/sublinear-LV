using DrWatson
@quickactivate
    
using DifferentialEquations
using Random, Distributions
using LinearAlgebra
using Plots


## generalized LV dynamics

function interaction_matrix(S, μ, σ, scaled = true, symm = false)
    if scaled
        A = rand(Normal(μ/S, σ/sqrt(S)), (S,S))
    else
        A = rand(Normal(μ, σ), (S,S))
    end
    if symm
        A = (A + A')/2
    end
    A[diagind(A)] .= 0
    return A
end

function ppart(x)
    return max(x, 1e-10)
end


function drift(u, p, t)
    return u.^p["k"] .- u.^2. / p["K"]^(2-p["k"]) .- u.*(p["A"]*u) .+ p["imm"]
end

function drift!(du, u, p, t)
        u .= ppart.(u)
        du .= u.^p["k"] .- u.^2. / p["K"]^(2-p["k"]) .- u.*(p["A"]*u) .+ p["imm"]
end

function jac!(J, u, p, t)
    J .= Diagonal(p["k"]*u.^(p["k"-1]) - 2u/p["K"]^(2-p["k"]) - p["A"]*u) - u.*p["A"] 
end

function diffusion!(du, u, p, t)
    u .= ppart.(u)
    du .= sqrt.(p["T"]*u)
end


## solving
blowup(threshold = 1e3) = DiscreteCallback((u, t, integrator) -> maximum(u) > threshold, terminate!)
collapse(fraction = 0.1) = DiscreteCallback((u, t, integrator) -> diversity(u) < fraction*length(u), terminate!)
converged(abstol = 1e-3) = TerminateSteadyState(abstol)

function  LV_problem(params; k = 1., x₀ = nothing, max_time = 1000.)
    
    p = copy(params)
    p["imm"] = fill(p["λ"], p["S"])
    p["A"] = interaction_matrix(p["S"], p["μ"], p["σ"], p["symm"], p["scaled"])
    
    if x₀ === nothing
        x₀ = rand(p["S"])
    end
    
    p["k"] = k

    if p["T"] == 0.
        pb = ODEProblem(drift!, x₀, (0., max_time), p);
    else
        pb = SDEProblem(drift!, diffusion!, x₀, (0., max_time), p);
    end
    
    return pb
end

function  LV_solver(params; k = 1., x₀ = nothing, termination = [blowup()], max_time = 1000.)
    p = copy(params)
    p["imm"] = fill(p["λ"], p["S"])
    p["A"] = interaction_matrix(p["S"], p["μ"], p["σ"], p["symm"], p["scaled"])
    if x₀ === nothing
        x₀ = rand(p["S"])
    end
    
    p["k"] = k

    if p["T"] == 0.
        pb = ODEProblem(drift!, x₀, (0., max_time), p);
    else
        pb = SDEProblem(drift!, diffusion!, x₀, (0., max_time), p);
    end

    p["sol"] = solve(pb, callback = CallbackSet(termination...))
    return p
end


function  LV_steady_state(params; k = 1., x₀ = nothing, termination = [blowup()], max_time = 1000.)
    p = copy(params)
    p["imm"] = fill(p["λ"], p["S"])
    p["A"] = interaction_matrix(p["S"], p["μ"], p["σ"], p["symm"], p["scaled"])
    if x₀ === nothing
        x₀ = rand(p["S"])
    end
    
    p["k"] = k

    if p["T"] == 0.
        pb = SteadyStateProblem(drift!, x₀, p);
    else
        error("Steady-state only possible at zero temperature")
    end

    solve(pb, callback = CallbackSet(termination...))
end

## plotting

function LV_compare_plot(params; max_time = 100)

    p = copy(params)
    p["imm"] = fill(p["λ"], p["S"])
    p["A"] = interaction_matrix(p["S"], p["μ"], p["σ"], p["symm"], p["scaled"])

    seed = rand(1:1000)
    Random.seed!(seed);
    IC = rand(p["S"])

    p["k"] = 1
    pb = SDEProblem(drift!, diffusion!, IC, (0., max_time), p);
    sol_lin = solve(pb);
    plin = plot(sol_lin, legend = false, title = "linear prod.");

    p["k"] = 0.75
    Random.seed!(seed);
    pb = SDEProblem(drift!, diffusion!, IC, (0., max_time), p);
    sol_sub = solve(pb);
    psub = plot(sol_sub, legend = false, title = "sublinear prod.");

    return plot(plin, psub)
end

## lyapunov spectrum

using DynamicalSystems
function  LV_lyap(params; k, x₀ = nothing, termination = [PositiveDomain(), blowup()], max_time = 1000.)
    p = copy(params)
    p["imm"] = fill(p["λ"], p["S"])
    p["A"] = interaction_matrix(p["S"], p["μ"], p["σ"], p["symm"], p["scaled"])
    if x₀ === nothing
        x₀ = rand(p["S"])
    end
    p["k"] = k

    if p["T"] == 0.
        pb = ContinuousDynamicalSystem(ODEProblem(drift!, x₀, (0., max_time), p));
    else
        error("Lyapunov spectrum only at zero T")
    end
    
    lyapunov(pb, max_time, callback = CallbackSet(termination...))
end