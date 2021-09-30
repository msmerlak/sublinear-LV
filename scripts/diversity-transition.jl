using DrWatson
@quickactivate

threshold = 1e3

include("../src/generalized-LV.jl")
using DataFrames, VegaLite, ThreadTools


function sim(p, k)
    sol = LV_solver(
        p, k = k, max_time = 1000., 
        termination = [blowup(threshold), converged()]
    )
    return (
        μ = p["μ"], 
        σ = p["σ"], 
        diversity = diversity(sol), 
        blowup = maximum(sol[end]) >= .9*threshold
        )
end

function sims(; S, T = 0, λ = 1e-3, symm = true, K)
    P = Dict{String, Any}("S" => S, "μ" => collect(-1:1:50), "σ" => collect(0:.05:5.), "λ" => λ, "T" => T, "K" => S, "symm" => symm, "scaled" => true, "det" => "yes")

    k = 1
    lin = tmap(p->sim(p, 1.), dict_list(P)) |> DataFrame
    tagsave(datadir("diversity-transition/", savename(@dict(S, K, k, T, λ), "csv")), lin)

    lin_plt = lin[lin.blowup .== false, :] |> @vlplot(:point, x=:μ, y=:σ, color=:diversity, title = "S = $S, k = $k, K = $K, T = $T") |> save(plotsdir("diversity-transition/", savename(@dict(S, K, k, T, λ), "csv")))

    k = 0.75
    sub = tmap(p->sim(p, k), dict_list(P)) |> DataFrame
    tagsave(datadir("diversity-transition/", savename(@dict(S, K, k, T, λ), "csv")), sub)
    sub_plt = sub[sub.blowup .== false, :] |> @vlplot(:point, x=:μ, y=:σ, color=:diversity, title = "S = $S, k = $k, K = $K, T = $T") |> save(plotsdir("diversity-transition/", savename(@dict(S, K, k, T, λ), "csv")))
end