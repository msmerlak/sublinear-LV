using DrWatson
@quickactivate

include("../src/generalized-LV.jl")
using DataFrames, VegaLite
using ThreadTools
## varying μ and σ at T = 0


function sim(p, k)
    sol = LV_solver(p, k = k, termination = [PositiveDomain(), blowup()], max_time = 200)
    return (success = sol.retcode, μ = p["μ"], σ = p["σ"], diversity = diversity(sol))
end


# symmetric
P = Dict{String, Any}("S" => 100, "μ" => collect(-1:.05:20), "σ" => collect(0:.05:.5), "λ" => 1e-2, "T" => .0, "K" => 100, "symm" => true, "scaled" => true, "det" => "yes")

lin = tmap(p->sim(p, 1.), dict_list(P)) |> DataFrame
sub = tmap(p->sim(p, 0.75), dict_list(P)) |> DataFrame


lin_plt = (lin[lin.success .== :Success, :] |> @vlplot(:point, x=:μ, y=:σ, color=:diversity, title = "k = 1")) |> save("phase-diag-lin.pdf")

sub[sub.success .== :Success, :] |> @vlplot(:point, x=:μ, y=:σ, color=:diversity, title = "k = 3/4") |> save("phase-diag-sublin.pdf")


# non-symmetric
P = Dict{String, Any}("S" => 100, "μ" => collect(-1:.05:20), "σ" => collect(0:.05:.5), "λ" => 1e-2, "T" => .0, "K" => 100, "symm" => false, "scaled" => true, "det" => "yes")

lin = tmap(p->sim(p, 1.), dict_list(P)) |> DataFrame
sub = tmap(p->sim(p, 0.75), dict_list(P)) |> DataFrame


lin_plt = (lin[lin.success .== :Success, :] |> @vlplot(:point, x=:μ, y=:σ, color=:diversity, title = "k = 1")) |> save("phase-diag-lin-nosymm.pdf")

sub[sub.success .== :Success, :] |> @vlplot(:point, x=:μ, y=:σ, color=:diversity, title = "k = 3/4") |> save("phase-diag-sublin-nosymm.pdf")

p = dict_list(P)[end]

p["μ"] = 1.
p["σ"] = .1
LV_solver(p, k = 1., termination = [PositiveDomain(), blowup()], max_time = 200.) |> plot
