using DrWatson
@quickactivate

using Distributed, ThreadTools
addprocs(Sys.CPU_THREADS-1)

@everywhere begin 
    threshold = 1e3
    include("../src/generalized-LV.jl")
    include("../src/metrics.jl")
    include("../src/simulations.jl")
end

using Plots

using DataFrames, VegaLite
using CSVFiles, ImageIO



for S in [10, 100, 500]
    println("S = $S")
    sims(S = S, K = 100, symm = true, scaled = false)
end

S = 100; T = 0.05; λ = 1e-2; symm = true; scaled = false
μ, σ = 50., 0.2

K = 100.

p = Dict{String, Any}("S" => S, "μ" => μ, "σ" => σ, "λ" => λ, "T" => T, "K" => K, "symm" => symm, "scaled" => scaled, "det" => "yes")

pp = LV_solver(
        p, k = 1., max_time = 2000., 
        termination = [blowup(threshold)]
    )

plot(pp["sol"], legend=false)



function F(u, p)
    return u.^p["k"] .- u.^2. / p["K"]^(2-p["k"]) .- u.*(p["A"]*u) .+ fill(p["λ"], p["S"])
end

f = x -> F(x, p)

using ForwardDiff
ForwardDiff.jacobian(x->F(x, pp), pp["sol"].u[end])