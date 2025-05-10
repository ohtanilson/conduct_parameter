using Distributed
Distributed.@everywhere include("../code/00setting_julia.jl")
Distributed.@everywhere include("../code/00functions.jl")
parameter = market_parameters_log()

estimation_methods = [
    #(:mpec,:theta_constraint, :slope_constraint, :equilibrium_constraint), 
    (:mpec,:no_constraint, :slope_constraint, :equilibrium_constraint),
    (:mpec,:theta_constraint, :non_constraint, :non_constraint),
    (:mpec,:no_constraint, :non_constraint, :non_constraint)
    ]
starting_value = :true_value
tol_level = :loose

## Estimate the parameters for each number of markets and the value of the standard deviation of the error terms

for estimation_method = estimation_methods
    for t = [100, 200, 1000, 1500], sigma =  [0.5, 1, 2]
        # Load the simulation data from the rds files
        filename_begin = "../conduct_parameter/output/data_loglinear_loglinear_n_"
        filename_end   = ".rds"

        if sigma == 1 || sigma == 2
            sigma = Int64(sigma)
        end
        filename = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_end

        data = load(filename)
        data = DataFrames.sort(data, [:group_id_k])
        # Set parameter values
        parameter = market_parameters_log(T = t, σ = sigma)
        
        # Estimation based on 2SLS
        @time estimation_result = iterate_estimation_nonlinear_2SLS(parameter, data, estimation_method, starting_value, tol_level)

        # Save the estimation result as csv file. The file is saved at "output" folder
        filename_estimation = "_"*String(estimation_method[1])*"_"*String(estimation_method[2])*"_"*String(estimation_method[3])

        filename_begin = "../conduct_parameter/output/parameter_hat_table_loglinear_loglinear_n_"
        filename_end   = ".csv"
        file_name = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_estimation*filename_end
        print("Simulate : $file_name \n")

        CSV.write(file_name, estimation_result[:,1:end-2], transform=(col, val) -> something(val, missing))

        histogram_distance_quantity = histogram((estimation_result[:,end]), bins = 100, title = "Distance between observed quantity and estimated quantity", ylabel = "Frequency", xlabel = "MSE", size = (1000, 600))

        display(histogram_distance_quantity)

        # Check if there are any non-missing values and if any of them are greater than 0
        if any(skipmissing(estimation_result[:,end-1]) .> 0)
            histogram_num_negative_Q = histogram(estimation_result[:,end-1], bins = 100, title = "Number of negative quantity", ylabel = "Frequency", xlabel = "Number of negative quantity", size = (1000, 600))
            display(histogram_num_negative_Q)
        end
    end

    println("\n")
    println("----------------------------------------------------------------------------------\n")
end


##
@everywhere function GMM_estimation_MPEC(T, Q, P, Z, Z_s, Z_d, X, X_s, X_d, α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3, θ_0, start_θ, start_γ,estimation_method::Tuple{Symbol, Symbol, Symbol, Symbol}, starting_value, tol_level)
    
    """
    Estimate the demand and supply parameter given a market simultaneously

    """
    if tol_level == :tight
        tol = 1e-15
        acceptable_tol = 1e-12
    elseif tol_level == :loose
        tol = 1e-6
        acceptable_tol = 1e-5
    end

    L = size(Z, 2)
    L_d = size(Z_d,2)
    L_s = size(Z_s,2)
    K_s = size(X_s, 2)
    K_d = size(X_d, 2)

    if estimation_method[1] == :mpec

        start_β = zeros(8)
        start_θ = 0
        if starting_value == :true_value
            start_β = [α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3]
            start_θ = θ_0
        elseif starting_value == :random
            start_β = [α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3] .+ rand(Uniform(-10, 10), 8)
            start_θ = θ_0 + rand(Uniform(-10, 1))
        end

        Ω = inv(Z' * Z)/T

        model = Model(Ipopt.Optimizer)
        set_optimizer_attribute(model, "tol", tol)
        set_optimizer_attribute(model, "max_iter", 1000)
        set_optimizer_attribute(model, "acceptable_tol", acceptable_tol)
        set_silent(model)
        @variable(model, β[k = 1:K_d+K_s-1], start = start_β[k])

        if estimation_method[2] == :theta_constraint
            @variable(model, 0 <= θ <= 1, start = start_θ)
        else
            @variable(model, θ, start = start_θ)
        end

        @variable(model, 0 <= MC[t = 1:T])
        for t = 1:T
            @NLconstraint(model, MC[t] == (1  - θ * (β[2] + β[3] * X[2*t,end]))* exp(P[t]))
        end

        ε = Any[];
        for t =1:T
            push!(ε, @NLexpression(model, P[t] - sum(β[k] * X[2*t-1,k] for k = 1:K_d) ))
            push!(ε, @NLexpression(model, log(MC[t]) - sum(β[k] * X[2*t,k] for k = K_d+1:K_d+K_s-1)))
        end
            
        if estimation_method[3] == :slope_constraint
            @constraint(model, β[K_d+2] >=0)                    # upward-sloping marginal cost
            for t = 1:T
                @NLconstraint(model, β[2] + β[3] * X_s[t, end] >= 0)                    # downward-sloping demand
            end
        end

        if estimation_method[4] == :equilibrium_constraint
            for t = 1:T
                @NLconstraint(model, 1 - θ * (β[2] + β[3] * X_s[t, end]) >= 0)          # Equilibrium constraint
            end
        end

        g = Any[];
        for l = 1:L
            push!(g, @NLexpression(model, sum(Z[t,l] * ε[t] for t = 1:2*T)))
        end
        @NLobjective(model, Min, sum( g[l] * Ω[l,k] * g[k] for l = 1:L, k = 1:L))
        optimize!(model)

        α_hat = value.(β)[1:K_d]
        γ_hat = value.(β)[K_d+1:end]
        θ_hat = value.(θ)
        ε_hat = value.(ε)

        num_negative_Q = 0
        dist = 0
        if termination_status(model) == MathOptInterface.LOCALLY_SOLVED
            num_negative_Q, dist = check_equilibrium_quantity(α_hat, γ_hat, θ_hat, ε_hat, Q, P, Z, Z_s, Z_d, X, X_s, X_d)
        else
            num_negative_Q = missing
            dist = missing
        end
        return α_hat, γ_hat, θ_hat, termination_status_code(termination_status(model)), num_negative_Q, dist
    end
end


@everywhere function check_equilibrium_quantity(α_hat, γ_hat, θ_hat, ε_hat, Q, P, Z, Z_s, Z_d, X, X_s, X_d)

    # Check if the reduced form quantity under the estimated parameter is consists of the data
    # for each t, if there is a reduced form that leads to negative quantity, check the number
    # Check also the distance between the observed data and reduced form quantity
    # the reduced form quantity is given by
    #log Q_t &= (\alpha_0 + \alpha_3 \log Y_t + \log (1 - \theta (\alpha_1 + \alpha_2 Z^{R}_{t})) - \gamma_0  -  \gamma_2 \log W_{t} - \gamma_3 \log R_t + \varepsilon^{d}_{t} - \varepsilon^{c}_{t})/(\gamma_1+ \alpha_1 + \alpha_2 Z^{R}_{t})

    # substitute the estimated parameter into the reduced form quantity

    # Devide the ε_hat into two parts
    # odd element is \varepsilon^{d}_{t}
    # even element is \varepsilon^{c}_{t}

    T = size(X_d, 1)

    ε_hat_d = ε_hat[mod.(1:2*T, 2) .== 1]
    ε_hat_c = ε_hat[mod.(1:2*T, 2) .== 0]

    log_Q_hat = (α_hat[1] .+ α_hat[4] .* X_d[:, end] .+ log.(1 .- θ_hat .* (α_hat[2] .+ α_hat[3] .* X_s[:, end])) .- γ_hat[1]  .-  γ_hat[3] .* X_s[:, 2] .- γ_hat[4] .* X_s[:, 3] .+ ε_hat_d[:] .- ε_hat_c[:])./(γ_hat[2] .+ α_hat[2] .+ α_hat[3] .* X_s[:,end])

    # check how many t are there that Q_hat is negative
    num_negative_Q = sum(exp.(log_Q_hat) .< 0)

    # check the distance between the observed data and reduced form quantity
    dist = sum((Q- log_Q_hat).^2)/T

    return num_negative_Q, dist
end