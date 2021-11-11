using DrWatson
@quickactivate

using Plots

threshold = 1e3
include("../src/generalized-LV.jl")

## multiple initial conditions
function  multi_eq(params; k = 1., x₀ = nothing, termination = [blowup(), converged()], max_time = 1000.)
    p = copy(params)
    p["imm"] = fill(p["λ"], p["S"])
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



S = 500; T = 0.0; λ = 0.; symm = true; scaled = true
μ, σ= 3., .5

k = .75

p = Dict{String, Any}("S" => S, "μ" => μ, "σ" => σ, "λ" => λ, "T" => T, "K" => K, "symm" => symm, "scaled" => scaled, "det" => "yes")

p["A"] = interaction_matrix(S, μ, σ, true, true)

p["K"] = 3.
plot(
        multi_eq(p, k = .75)["sol"],
        legend = false

    )

