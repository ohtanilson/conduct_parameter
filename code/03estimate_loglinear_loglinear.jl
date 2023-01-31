#---------------------------------------------------------------------------------------------------------

function termination_status_code(status)

    # See https://jump.dev/JuMP.jl/stable/moi/reference/models/#MathOptInterface.TerminationStatusCode
    
    status_code = [
        JuMP.MathOptInterface.OPTIMAL,                   # 1
        JuMP.MathOptInterface.INFEASIBLE,                # 2
        JuMP.MathOptInterface.DUAL_INFEASIBLE,           # 3
        JuMP.MathOptInterface.LOCALLY_SOLVED,            # 4
        JuMP.MathOptInterface.LOCALLY_INFEASIBLE,        # 5
        JuMP.MathOptInterface.INFEASIBLE_OR_UNBOUNDED,   # 6
        JuMP.MathOptInterface.ALMOST_OPTIMAL,            # 7
        JuMP.MathOptInterface.ALMOST_INFEASIBLE,         # 8
        JuMP.MathOptInterface.ALMOST_DUAL_INFEASIBLE,    # 9
        JuMP.MathOptInterface.ALMOST_LOCALLY_SOLVED,     # 10
        JuMP.MathOptInterface.ITERATION_LIMIT,           # 11
        JuMP.MathOptInterface.TIME_LIMIT,                # 12
        JuMP.MathOptInterface.NODE_LIMIT,                # 13
        JuMP.MathOptInterface.SOLUTION_LIMIT,            # 14
        JuMP.MathOptInterface.MEMORY_LIMIT,              # 15
        JuMP.MathOptInterface.OBJECTIVE_LIMIT,           # 16
        JuMP.MathOptInterface.NORM_LIMIT,                # 17
        JuMP.MathOptInterface.OTHER_LIMIT,               # 18
        JuMP.MathOptInterface.SLOW_PROGRESS,             # 19
        JuMP.MathOptInterface.NUMERICAL_ERROR,           # 20
        JuMP.MathOptInterface.INVALID_MODEL,             # 21
        JuMP.MathOptInterface.INVALID_OPTION,            # 22
        JuMP.MathOptInterface.INTERRUPTED,               # 23
        JuMP.MathOptInterface.OTHER_ERROR                # 24
    ]

    for i = eachindex(status_code)
        if (status == status_code[i]) == true
            return i
            break
        end
    end
end

#---------------------------------------------------------------------------------------------

function GMM_estimation_separate(T, Q, P, Z, Z_s, Z_d, X, X_s, X_d, parameter, estimation_method::Tuple{Symbol, Symbol, Symbol})
    
    """
    Estimate the demand and supply parameter given a market
    The demand parameters are destimated by the OLS.
    The supply parameters and the conduct parameter are estimated by the GMM.
    """

    @unpack θ, start_θ, start_γ = parameter

    QZ_hat = Z_d * inv(Z_d' * Z_d) * Z_d' * (Z_d[:,2] .* Q)
    Q_hat = Z_d * inv(Z_d' * Z_d) * Z_d' *  Q
    X_dd = hcat(ones(T), -Q_hat, -QZ_hat, X_d[:,end])
    α_hat = inv(X_dd' * X_dd) * (X_dd' * P)

    L = size(Z, 2)
    L_d = size(Z_d,2)
    L_s = size(Z_s,2)
    K_s = size(X_s, 2)
    K_d = size(X_d, 2)

    sample_violatiton_index = Int64[]

    # Pick up the index of market under which the inside of the log has a negative value
    for t = 1:T
        if 1 .- θ .*(α_hat[2] .+ α_hat[3] .* X_s[t,end]) <= 0
            push!(sample_violatiton_index, t)
        end
    end

    # If all markets do not satisfy the assumption, stop the supply estimation and rerutns missing values
    if length(sample_violatiton_index) == T

        γ_hat = repeat([missing], K_s-1)
        θ_hat = missing

        return α_hat, γ_hat, θ_hat, missing

    else
        
        # Drop the samples that violate the assumption
        if  1 <= length(sample_violatiton_index)

            sample_index = setdiff([1:T;], sample_violatiton_index)

            Z_s = Z_s[sample_index, :]
            X_s = X_s[sample_index, :]
            T = length(sample_index)
            L_s = size(Z_s,2)
            K_s = size(X_s, 2)
        end

        #Check if the weight matrix can be obtained
        if rank(Z_s' * Z_s) == L_s
            
            Ω = inv(Z_s' * Z_s)/T

            model = JuMP.Model(Ipopt.Optimizer)
            JuMP.set_optimizer_attribute(model, "tol", 1e-15)
            JuMP.set_optimizer_attribute(model, "max_iter", 1000)
            JuMP.set_optimizer_attribute(model, "acceptable_tol", 1e-12)
            JuMP.set_silent(model)
            JuMP.@variable(model, γ[k = 1:K_s-1], start = start_γ[k])

            if estimation_method[3] == :theta_constraint
                JuMP.@variable(model, 0 <= θ <= 1, start = start_θ)
            else
                JuMP.@variable(model, θ, start = start_θ)
            end

            r = Any[];
            g = Any[];
            for t =1:T
                push!(r, JuMP.@NLexpression(model, P[t]- sum(γ[k] * X_s[t,k] for k = 1:K_s-1) + log(1 - θ *(α_hat[2] + α_hat[3] * X_s[t, end]) ) ) )
            end

            for l = 1:L_s
                push!(g, JuMP.@NLexpression(model, sum(Z_s[t,l] * r[t] for t = 1:T)))
            end

            if estimation_method[2] == :log_constraint

                for t = 1:T
                    JuMP.@NLconstraint(model, 0 <= 1 - θ *(α_hat[2] + α_hat[3] * X_s[t, end]) )
                end
            end

            JuMP.@NLobjective(model, Min, sum( g[l] *Ω[l,k] * g[k] for l = 1:L_s, k = 1:L_s))

            JuMP.optimize!(model)
            
            γ_hat = JuMP.value.(γ)
            θ_hat = JuMP.value.(θ)
        else

            # If the weight matrix can not be obtained, return missing values
            γ_hat = repeat([missing], K_s-1)
            θ_hat = missing

            return α_hat, γ_hat, θ_hat, missing
        end

        # Check if the supply estimation result satisfies the assumption
        if  sum(1 .- θ_hat .*(α_hat[2] .+ α_hat[3] .* X_s[:,end]) .<= 0) == 0

            return α_hat, γ_hat, θ_hat, termination_status_code(JuMP.termination_status(model))
        else 
            error("The estimation result violates the model assumption ")
        end
    end

end

function estimation_nonlinear_2SLS(parameter, data, estimation_method::Tuple{Symbol, Symbol, Symbol})
    """
    Given data, reshape the data and pass it to the function that implement the GMM estimation

    """

    @unpack T = parameter

    Q  = data.logQ
    w  = data.w
    r  = data.r
    z  = data.z
    iv_w = data.iv_w
    iv_r = data.iv_r
    p  = data.logP
    y  = data.y

    iv = hcat(iv_w, iv_r)

    X_d = []
    X_s = []
    X   = []
    Z   = []
    Z_d = []
    Z_s = []
    P   = []

    for t = 1:T

        Z_dt = vcat(1, z[t], iv[t,:], y[t])
        Z_st = vcat(1, z[t], log(w[t]), log(r[t]), y[t])
        Z_t = [Z_dt zeros(length(Z_dt));  zeros(length(Z_st)) Z_st]'

        X_dt = vcat(1, -Q[t], -z[t].*Q[t], y[t])
        X_st = vcat(1, Q[t], log(w[t]), log(r[t]), z[t])
        X_t  = [X_dt zeros(length(X_dt));  zeros(length(X_st)) X_st]'

        push!(P, p[t])
        push!(X_d, X_dt')
        push!(X_s, X_st')
        push!(X, X_t)
        push!(Z, Z_t)
        push!(Z_d, Z_dt')
        push!(Z_s, Z_st')
    end

    Z   = reduce(vcat,(Z))
    X   = reduce(vcat,(X))
    X_d = reduce(vcat,(X_d))
    X_s = reduce(vcat,(X_s))
    Z_d = reduce(vcat,(Z_d))
    Z_s = reduce(vcat,(Z_s))

    if estimation_method[1] == :separate
        α_hat, γ_hat, θ_hat, status = GMM_estimation_separate(T, Q, P, Z, Z_s, Z_d, X, X_s, X_d, parameter, estimation_method)
    else
        α_hat, γ_hat, θ_hat, status = GMM_estimation_simultaneous(T, Q, P, Z, Z_s, Z_d, X, X_s, X_d, parameter, estimation_method)
    end

    return α_hat, γ_hat, θ_hat, status
end


function simulation_nonlinear_2SLS(parameter, data, estimation_method::Tuple{Symbol, Symbol, Symbol})

    """
    Given the simulation data, run the estimation in each simulation index s = 1,..., 1000, and store the simulation results as a DataFrame file.

    """

    @unpack T, S, σ = parameter 

    α_est  = Vector{Float64}[]
    γ_est  = Vector{Union{Missing, Float64}}[]
    θ_est  = Union{Missing, Float64}[]
    status = Union{Missing, Int64}[]

    for s = 1:S
        data_s = data[(s-1)*T+1:s*T,:]
        α_est_s, γ_est_s, θ_est_s, status_s = estimation_nonlinear_2SLS(parameter, data_s, estimation_method)
    
        push!(α_est, α_est_s)
        push!(γ_est, γ_est_s)
        push!(θ_est, θ_est_s)
        push!(status,status_s)
    end
    
    α_est = reduce(vcat, α_est')
    γ_est = reduce(vcat, γ_est')
    θ_est = reduce(vcat, θ_est')
    status = reduce(vcat, status)

    status_indicator = []

    for i = eachindex(status)
        if status[i] !== missing
            if status[i] <= 7
                push!(status_indicator, 1)
            else
                push!(status_indicator, 0)
            end
        else
            push!(status_indicator, missing)
        end
    end

    estimation_result = DataFrame(
    T = T,
    σ = σ,
    α_0 = α_est[:,1],
    α_1 = α_est[:,2],
    α_2 = α_est[:,3],
    α_3 = α_est[:,4],
    γ_0 = γ_est[:,1],
    γ_1 = γ_est[:,2],
    γ_2 = γ_est[:,3],
    γ_3 = γ_est[:,4],
    θ = θ_est,
    status = status,
    status_indicator = status_indicator)

    return estimation_result
end


#--------------------------------------------------------------------------------------------------------------

# Estimate the parameters for each number of markets and the value of the standard deviation of the error terms

# Estimation start from with default starting values
for estimation_method = estimation_methods
    for t = [50, 100, 200, 1000], sigma =  [0.001, 0.5, 1, 2]

        # Load the simulation data from the rds files
        filename_begin = "../conduct_parameter/output/data_loglinear_loglinear_n_"
        filename_end   = ".rds"

        if sigma == 1 || sigma == 2
            sigma = Int64(sigma)
        end

        filename = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_end

        data = load(filename)
        data = sort(data, [:group_id_k])

        #Uncomment the following lines to load the simulation data from the csv files
            #filename_begin = "../conduct_parameter/output/data_loglinear_loglinear_n_"
            #filename_end   = ".csv"
            #file_name = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_end
            #data = DataFrame(CSV.File(file_name))

        # Set parameter values
        parameter = market_parameters_log(T = t, σ = sigma)
        
        # Estimation based on 2SLS
        @time estimation_result = simulation_nonlinear_2SLS(parameter, data, estimation_method)

        # Save the estimation result as csv file. The file is saved at "output" folder
        filename_estimation = "_"*String(estimation_method[1])*"_"*String(estimation_method[2])*"_"*String(estimation_method[3])

        filename_begin = "../conduct_parameter/output/parameter_hat_table_loglinear_loglinear_n_"
        filename_end   = ".csv"
        file_name = filename_begin*string(t)*"_sigma_"*string(sigma)*filename_estimation*filename_end

        CSV.write(file_name, estimation_result, transform=(col, val) -> something(val, missing))

    end
    println("\n")
    println("----------------------------------------------------------------------------------\n")
end