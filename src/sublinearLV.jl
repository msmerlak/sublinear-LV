using DifferentialEquations
using Random, Distributions
using Plots; plotly
using LinearAlgebra

function interaction_matrix(; S, μ, σ, symm = false)
    A = rand(Normal(μ/S, σ/sqrt(S)), (S,S))
    if symm
        A = (A + A')/2
    end
    A[diagind(A)] .= 0
    return A
end

function ppart(x)
    return max(x, 0)
end


function f!(du, u, p, t)
    v = ppart.(u)
    du .= v.*(v.^(p[:k]-1) - (v/p[:K]).^2 - p[:A]*v) + p[:imm]
end

function g!(du, u, p, t)
    du .= sqrt.(p[:T]*ppart.(u))
end

IC = rand(15)

function LV_plot(params; max_time = 100)

    p = copy(params)
    p[:imm] = fill(p[:λ], p[:S])
    p[:A] = interaction_matrix(S = p[:S], μ = p[:μ], σ = p[:σ], symm = p[:symm])

    seed = rand(1:1000)
    Random.seed!(seed);
    IC = rand(p[:S])

    p[:k] = 1
    pb = SDEProblem(f!, g!, IC, (0., max_time), p);
    sol = solve(pb);
    plin = plot(sol, legend = false, title = "linear prod.");

    p[:k] = 0.75
    Random.seed!(seed);
    pb = SDEProblem(f!, g!, IC, (0., max_time), p);
    sol = solve(pb);
    psub = plot(sol, legend = false, title = "sublinear prod.");

    plot(plin, psub)
end
