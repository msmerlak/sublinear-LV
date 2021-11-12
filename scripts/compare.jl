using DrWatson
@quickactivate

include(srcdir("plotting.jl"))

p = Dict{String, Any}(
    "S" => 100,
    "μ" => 1.,
    "σ" => .05,
    "symm" => false,
    "scaled" => true,
    "λ" => 0.,
    "K" => 10.,
    "T" => 0., 
    "ξ" => .5,
)


plot_sol(p)


savefig(current(), plotsdir("2021-11-11/", "asymmetric.png"))

s = sim2(p)

using StatsBase
lin_cdf = ecdf(s.lin.abundances)
sub_cdf = ecdf(s.sub.abundances)


plot(x->1-lin_cdf(x), xlims = (1e-4, maximum(s.lin.abundances)), xaxis = :log, yaxis = :log)
plot(x->1-sub_cdf(x), xlims = (1e-4, maximum(s.sub.abundances)), xaxis = :log)

plot(s.lin.trajectory)

histogram(log10.(s.lin.abundances[s.lin.abundances .> 1e-3]))
histogram(log10.(s.sub.abundances))




plts = plot_K(P)
savefig(current(), plotsdir("increase-K.png"))