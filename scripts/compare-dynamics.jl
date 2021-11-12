using DrWatson
@quickactivate
include("../src/generalized-LV.jl")

using ProgressMeter

## Phase transition along the μ axis

p = Dict{Symbol, Any}(:S => 100, :μ => .2, :σ => .1, :λ => 0.01, :T => .001, :K => 100, :symm => false, :scaled => true)

#P = Plots.Plot{Plots.GRBackend}[]
p[:μ] = 0.
plt = LV_plot(p, max_time = 200.)

for μ=0:1:10
    p[:μ] = μ
    println(μ)
    plt = LV_compare_plot(p, max_time = 100.)
    plt = plot(plt, plot_title = savename(p), titlefont=font(10,"Computer Modern"));
#    push!(P, plt)
    png(plt, savename("plots/", p))
end