using Distributed
Distributed.@everywhere include("../code/00setting_julia.jl")
Distributed.@everywhere include("../code/00functions.jl")
parameter = market_parameters_log()

estimation_methods = [
    (:mpec,:theta_constraint, :slope_constraint, :equilibrium_constraint), 
    (:mpec,:no_constraint, :slope_constraint, :equilibrium_constraint),
    (:mpec,:theta_constraint, :non_constraint, :non_constraint),
    (:mpec,:no_constraint, :non_constraint, :non_constraint)
    ]
starting_value = :true_value
tol_level = :loose

## Estimate the parameters for each number of markets and the value of the standard deviation of the error terms

for estimation_method = estimation_methods
    for t = [100], sigma =  [0.5]
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

        CSV.write(file_name, estimation_result, transform=(col, val) -> something(val, missing))
    end
    println("\n")
    println("----------------------------------------------------------------------------------\n")
end



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
            @NLconstraint(model, exp(P[t]) == MC[t] + θ * (β[2] + β[3] * X[2*t,end])* exp(P[t]))
            #@NLconstraint(model, (exp(P[t]) - MC[t])/((β[2] + β[3] * X[2*t,end])* exp(P[t])) == θ)
        end

        r = Any[];
        for t =1:T
            push!(r, @NLexpression(model, P[t] - sum(β[k] * X[2*t-1,k] for k = 1:K_d) ))
            push!(r, @NLexpression(model, log(MC[t]) - sum(β[k] * X[2*t,k] for k = K_d+1:K_d+K_s-1)))
        end
            
        if estimation_method[3] == :slope_constraint
            #@constraint(model, β[K_d+2] >=0)                    # upward-sloping marginal cost
            for t = 1:T
                @NLconstraint(model, β[2] + β[3] * X_s[t, end] + β[K_d+2] >= 0)                    # downward-sloping demand
            end
        end

        if estimation_method[4] == :equilibrium_constraint
            for t = 1:T
                @NLconstraint(model, 1 - θ * (β[2] + β[3] * X_s[t, end]) >= 0)          # Equilibrium constraint
            end
        end

        g = Any[];
        for l = 1:L
            push!(g, @NLexpression(model, sum(Z[t,l] * r[t] for t = 1:2*T)))
        end
        @NLobjective(model, Min, sum( g[l] * Ω[l,k] * g[k] for l = 1:L, k = 1:L))
        optimize!(model)

        α_hat = value.(β)[1:K_d]
        γ_hat = value.(β)[K_d+1:end]
        θ_hat = value.(θ)

        return α_hat, γ_hat, θ_hat, termination_status_code(termination_status(model))

    end
end
