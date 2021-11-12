
using DrWatson
@quickactivate

include(srcdir("generalized-LV.jl"))
include(srcdir("simulations.jl"))
using Plots

## plotting

function plot_params(sol)
    A = sol.params["A"]
    
    int_matrix = heatmap(A, title = "μ = $(params["μ"]), σ = $(params["σ"])")
    
    surv = sol.abundances .> 0
    surv_matrix = heatmap(A[surv, surv], title = "μ = $(params["μ"]), σ = $(params["σ"])")

    grow_rates = histogram(sol.params["r"], xlabel = "growth rate rᵢ", legend = false)
    plot(int_matrix, surv_matrix)
end

function plot_sol(params)

    lin, sub = sim2(params)
    plin = plot(lin.trajectory, legend = false, title = "diversity = $(round(lin.diversity))", color = :skyblue1, alpha = 0.9);
    psub = plot(sub.trajectory, legend = false, title = "diversity = $(round(sub.diversity))", color = :orange, alpha = 0.9);

    pspec = scatter(lin.spectrum, color = :skyblue1, alpha = 0.5, aspect_ratio = 1, label = "k = 1", xlabel = "Re(eigenvalue)", ylabel = "Im(eigenvalue)")
    scatter!(pspec, sub.spectrum, color = :orange, alpha = 0.5, aspect_ratio = 1, label = "k = 3/4")

    lin_ab = lin.abundances[lin.abundances .> 0]
    hist = histogram(log10.(lin_ab), color = :skyblue1, alpha = 0.5, xlabel = "log10(abundance)", ylabel = "count", label = "k = 1")

    sub_ab = sub.abundances[sub.abundances .> 0]
    histogram!(log10.(sub_ab), color = :orange, alpha = 0.5, label = "k = 3/4")

    sc = scatter(lin.abundances, sub.abundances)
    return plot(plin, psub, pspec, hist)
end


function plot_K(params)
    plots = []
    for K in params["K"]
        p = deepcopy(params)
        p["K"] = K
        lin, sub = sim2(p)
        plin = plot(lin.trajectory, legend = false, ylims = (0, 1.1*K), size = (100, 40), xlabel = "")
        hline!([K], color=:black)
        psub = plot(sub.trajectory, legend = false, size = (1000, 400), xlabel = "")
        hline!([K], color=:black)
        push!(plots, plot(plin, psub))
    end
    #return plots
    plot(plots..., layout = (length(plots), 1))
end

