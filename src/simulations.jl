using DrWatson
@quickactivate

include(srcdir("generalized-LV.jl"))
include(srcdir("metrics.jl"))

using ForwardDiff

function sim(params, k; lyapunov = true, max_time = 1000, termination = [blowup(), converged()])
    sol, params = LV_solver(
        params, k = k, max_time = max_time, termination = termination
    )

    params["k"] = k
    jac = ForwardDiff.jacobian(u -> F(u, params), sol[end])

    if lyapunov
        λ = LV_lyap(
            params, k = k, max_time = max_time, termination = termination
        )
    end
    return (
        params = params,
        status = sol.retcode ,
        diversity = diversity(sol), 
        richness = richness(sol), 
        abundances = sol[end],
        trajectory = sol,
        jacobian = jac,
        spectrum = eigen(jac).values,
        largest_lyapunov_exponent = lyapunov ? λ : nothing
        )
end

function sim2(params)
    seed = rand(UInt)
    
    Random.seed!(seed)
    lin = sim(params, 1)

    Random.seed!(seed)
    sub = sim(params, .75)

    return (lin = lin, sub = sub)
end
