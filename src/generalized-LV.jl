using DifferentialEquations
using Random, Distributions
using LinearAlgebra



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
    return max(x, 0)
end


function drift!(du, u, p, t)
        v = ppart.(u)
        du .= v.^p["k"] .- v.^2. / p["K"]^(2-p["k"]) .- v.*(p["A"]*v) .+ p["imm"]
end

function diffusion!(du, u, p, t)
    du .= sqrt.(p["T"]*ppart.(u))
end

## diversity

function diversity(u::Vector)
    v = u./norm(u, 1)
    return exp(-sum(v.*log.(v)))
end

diversity(sol::RODESolution) = diversity(sol.u[end])

## solving
blowup() = DiscreteCallback((u, t, integrator) -> maximum(u) > 1e4, terminate!)
collapse() = DiscreteCallback((u, t, integrator) -> diversity(u) < length(u)/10, terminate!)

function  LV_solver(params; k, x₀ = nothing, termination = [nothing], max_time = 1000.)
    p = copy(params)
    p["imm"] = fill(p["λ"], p["S"])
    p["A"] = interaction_matrix(p["S"], p["μ"], p["σ"], p["symm"], p["scaled"])
    if x₀ == nothing
        x₀ = rand(p["S"])
    end
    p["k"] = k
    pb = SDEProblem(drift!, diffusion!, x₀, (0., max_time), p);

    cb = CallbackSet(termination...)

    solve(pb, callback = cb)
end



## plotting

function LV_compare_plot(params; max_time = 100)

    p = copy(params)
    p[:imm] = fill(p[:λ], p[:S])
    p[:A] = interaction_matrix(S = p[:S], μ = p[:μ], σ = p[:σ], symm = p[:symm], scaled = p[:scaled])

    seed = rand(1:1000)
    Random.seed!(seed);
    IC = rand(p[:S])

    p[:k] = 1
    pb = SDEProblem(drift!, diffusion!, IC, (0., max_time), p);
    sol_lin = solve(pb);
    plin = plot(sol_lin, legend = false, title = "linear prod.");

    p[:k] = 0.75
    Random.seed!(seed);
    pb = SDEProblem(drift!, diffusion!, IC, (0., max_time), p);
    sol_sub = solve(pb);
    psub = plot(sol_sub, legend = false, title = "sublinear prod.");

    return plot(plin, psub)
end
