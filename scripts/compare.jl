using DrWatson
@quickactivate

include(srcdir("simulations.jl"))
include(srcdir("plotting.jl"))


using Dates
day = string(today())


p = Dict{String, Any}(
    "S" => 100,
    "μ" => 10.,
    "σ" => .5,
    "symm" => true,
    "scaled" => true,
    "λ" => 0.,
    "K" => 100.,
    "T" => 0., 
    "ξ" => .0,
)

sols = sim2(p);
plot_sol(sols)

savefig(plotsdir(day * "/example-plot.pdf"))

p = Dict{String, Any}(
    "S" => 50,
    "μ" => 0.1,
    "σ" => 0.01,
    "symm" => false,
    "scaled" => true,
    "λ" => 0.,
    "K" => 20.,
    "T" => 0., 
    "ξ" => .5,
)

@time sols = sim2(p);
plot_sol(sols)