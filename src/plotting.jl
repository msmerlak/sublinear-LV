
using DrWatson
@quickactivate


using Plots, PrettyTables

## plotting

function plot_params(sol)

    A = sol.params["A"]
    
    int_matrix = heatmap(A, xlabel = "interaction coefficients")
    
    #surv = sol.abundances .> 1e-3
    #surv_matrix = heatmap(A[surv, surv], title = "μ = $(sol.params["μ"]), σ = $(sol.params["σ"])")

    grow_rates = histogram(sol.params["r"], xlabel = "growth rates", legend = false)
    plot(int_matrix, grow_rates)
end

function print_metrics(sols)
    lin, sub = sols
    metrics = [
        "diversity" lin.diversity   sub.diversity;
        "dominant eigenvalue" argmax(real, lin.spectrum) argmax(real, sub.spectrum);
        "LLE" lin.largest_lyapunov_exponent sub.largest_lyapunov_exponent
    ]
    return pretty_table(metrics, header = ["", "k = 1", "k = 3/4"])

end

function plot_sol(sols)

    print_metrics(sols)

    lin, sub = sols

    p1 = plot_params(lin)
    
    plin = plot(lin.trajectory, legend = false, color = :skyblue1, alpha = 0.9, xlabel = "time series (linear prod)");
    hline!([lin.params["K"]], color=:black)
    psub = plot(sub.trajectory, legend = false, color = :orange, alpha = 0.9, xlabel = "time series (sublinear prod)");
    hline!([lin.params["K"]], color=:black)
    p2 = plot(plin, psub)



    hist = histogram(lin.abundances, color = :skyblue1, alpha = 0.5, legend = false, xlabel = "abundance", label = "k = 1")
    histogram!(hist, sub.abundances, color = :orange, alpha = 0.5, legend = false, label = "k = 3/4")

    lin_ab = lin.abundances[lin.abundances .> 0]
    loghist = histogram(log10.(lin_ab), color = :skyblue1, alpha = 0.5, legend = false, xlabel = "log10(abundance)", label = "k = 1")

    sub_ab = sub.abundances[sub.abundances .> 0]
    histogram!(loghist, log10.(sub_ab), color = :orange, alpha = 0.5, legend = false, label = "k = 3/4")
    
    p3 = plot(hist, loghist)


    spectrum = scatter(lin.spectrum, color = :skyblue1, alpha = 0.5, size = (1000, 1000), xlabel = "Jacobian spectrum", ylabel = "", legend = false)
    scatter!(spectrum, sub.spectrum, color = :orange, alpha = 0.5, size = (1000, 1000), legend = false, xlabel = "Jacobian spectrum", ylabel = "")
    vline!([0], color = :black)

    diversity = bar([lin.diversity, sub.diversity], xlabel = "diversity", legend = false, alpha = 0.5, color = [:skyblue1, :orange])
    LLE = bar([- lin.largest_lyapunov_exponent, - sub.largest_lyapunov_exponent], xlabel = "- LLE", legend = false, alpha = 0.5, color = [:skyblue1, :orange])
    barplots = plot(diversity, LLE)

    p4 = plot(spectrum, barplots)
    # sc = scatter(lin.abundances, sub.abundances)
    return plot(p1, p2, p3, p4, layout = (4, 1))
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

