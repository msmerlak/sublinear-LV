include("../src/sublinearLV.jl")
using DrWatson

p = Dict{Symbol, Any}(:S => 100, :μ => 1, :σ => .1, :λ => 0.1, :T => .1, :K => 50, :symm => true)


plt = LV_plot(p, max_time = 100)
plt = plot(plt, plot_title = savename(p), titlefont=font(10,"Computer Modern"))
png(plt, savename("plots/", p))

p = Dict{Symbol, Any}(:S => 100, :μ => 1, :σ => .1, :λ => 0.1, :T => .01, :K => 50, :symm => false)

plt = LV_plot(p, max_time = 500)
plt = plot(plt, plot_title = savename(p), titlefont=font(10,"Computer Modern"))
png(plt, savename("plots/", p))
