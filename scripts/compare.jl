using DrWatson
@quickactivate

include(srcdir("simulations.jl"))
include(srcdir("plotting.jl"))


using Dates
day = string(today())


p = Dict{String, Any}(
    "S" => 100,
    "μ" => 1.,
    "σ" => .1,
    "symm" => false,
    "scaled" => true,
    "λ" => 0.,
    "K" => 20.,
    "T" => 0., 
    "ξ" => .5,
)

sols = sim2(p)
plot_sol(sols)

savefig(plotsdir(day * "/example-plot.pdf"))