function sim(p, k)
    p = LV_solver(
        p, k = k, max_time = 1000., 
        termination = [blowup(threshold), converged()]
    )
    return (
        μ = p["μ"],
        σ = p["σ"], 
        diversity = diversity(p["sol"]), 
        richness = richness(p["sol"]), 
        blowup = maximum(p["sol"][end]) >= .9*threshold
        )
end



# function sims(; S, T = 0, λ = 1e-3, symm = true, K, scaled = true)
#     P = Dict{String, Any}("S" => S, "μ" => collect(-1:1:50), "σ" => collect(0:.05:5.), "λ" => λ, "T" => T, "K" => S, "symm" => symm, "scaled" => scaled, "det" => "yes")

#     k = 1
#     lin = pmap(p->sim(p, 1.), dict_list(P)) |> DataFrame
#     save(datadir("diversity-transition/", savename(@dict(S, K, k, T, λ, symm, scaled), "csv")), lin)

#     lin[lin.blowup .== false, :] |> @vlplot(:point, x=:μ, y=:σ, xaxis = :log, yaxis = :log, color=:diversity, title = "S = $S, k = $k, K = $K, T = $T, symm = $symm, scaled = $scaled") |> save(plotsdir("diversity-transition/", savename(@dict(S, K, k, T, λ, symm, scaled), "png")))

#     k = 0.75
#     sub = pmap(p->sim(p, k), dict_list(P)) |> DataFrame
#     save(datadir("diversity-transition/", savename(@dict(symm, scaled), "csv")), sub)
#     sub_plt = sub[sub.blowup .== false, :] |> @vlplot(:point, x=:μ, y=:σ,xaxis = :log, yaxis = :log, color=:diversity, title = "S = $S, k = $k, K = $K, T = $T, symm = $symm, scaled = $scaled") |> save(plotsdir("diversity-transition/", savename(@dict(S, K, k, T, λ, symm, scaled), "png")))
# end