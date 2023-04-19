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

function GMM_estimation_separate(T, Q, P, Z, Z_s, Z_d, X, X_s, X_d, parameter, estimation_method::Tuple{Symbol, Symbol, Symbol}, start_value, tol_level)

    """
    Estimate the demand and supply parameter given a market
    The demand parameters are destimated by IV.
    The supply parameters and the conduct parameter are estimated by the GMM.
    """

    @unpack γ_0, γ_1, γ_2, γ_3, θ, start_θ, start_γ = parameter
    # first stage
    QZ_hat = Z_d * inv(Z_d' * Z_d) * Z_d' * (Z_d[:,2] .* Q)
    Q_hat = Z_d * inv(Z_d' * Z_d) * Z_d' *  Q
    # second stage
    X_dd = hcat(ones(T), -Q_hat, -QZ_hat, X_d[:,end])
    α_hat = inv(X_dd' * X_dd) * (X_dd' * P)

    L = size(Z, 2)
    L_d = size(Z_d,2)
    L_s = size(Z_s,2)
    K_s = size(X_s, 2)
    K_d = size(X_d, 2)

    sample_violatiton_index = Int64[]

    if tol_level == :tight
        tol = 1e-15
        acceptable_tol = 1e-12
    elseif tol_level == :loose
        tol = 1e-6
        acceptable_tol = 1e-5
    end

    if start_value == :true
        start_θ = θ
        start_γ = [γ_0, γ_1, γ_2, γ_3]

    elseif start_value == :random
        start_γ = [γ_0, γ_1, γ_2, γ_3] .+ rand(Uniform(-20, 20), 4)
        start_θ = 10
        while sum(1 .- start_θ .*(α_hat[2] .+ α_hat[3] .* X_s[:,end]) .<= 0) != 0
            start_θ = θ + rand(Uniform(-10, 1))
        end
    end

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
        if LinearAlgebra.rank(Z_s' * Z_s) == L_s

            Ω = inv(Z_s' * Z_s)/T

            model = JuMP.Model(Ipopt.Optimizer)
            JuMP.set_optimizer_attribute(model, "tol", tol)
            JuMP.set_optimizer_attribute(model, "max_iter", 1000)
            JuMP.set_optimizer_attribute(model, "acceptable_tol", acceptable_tol)
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

            # If the weight matrix cannot be obtained, return missing values
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


function GMM_estimation_simultaneous(T, Q, P, Z, Z_s, Z_d, X, X_s, X_d, parameter, estimation_method::Tuple{Symbol, Symbol, Symbol}, start_value, tol_level)
    
    """ 
    Estimate the demand and supply parameter given a market simultaneously
    The 
    
    """
    @unpack γ_0, γ_1, γ_2, γ_3, θ, start_θ, start_γ = parameter

    L = size(Z, 2)
    K_s = size(X_s, 2)
    K_d = size(X_d, 2)

    # The wight function for the GMM estimation
    Ω = inv(Z' * Z)/T

    if tol_level == :tight
        tol = 1e-15
        acceptable_tol = 1e-12
    elseif tol_level == :loose
        tol = 1e-6
        acceptable_tol = 1e-5
    end

    if start_value == :true
        start_β = [α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3]
        start_θ = θ

    elseif start_value == :random
        start_β = [α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3] .+ rand(Uniform(-10, 10), 8)
        while sum(1 .- start_θ .*(α_1[2] .+ α_2 .* X[:,end]) .<= 0) != 0
            start_θ = θ + rand(Uniform(-10, 1))
        end
    end

    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "tol", tol)
    set_optimizer_attribute(model, "max_iter", 1000)
    set_optimizer_attribute(model, "acceptable_tol", acceptable_tol)
    set_silent(model)
    @variable(model, β[k = 1:K_d+K_s-1])

    if estimation_method[3] == :theta_constraint
        @variable(model, 0 <= θ <= 1)
    else
        @variable(model, θ)
    end

    r = Any[];
    g = Any[];
    for t =1:T
        push!(r, @NLexpression(model, P[t]- sum(β[k] * X[2*t-1,k] for k = 1:K_d) ))
        push!(r, @NLexpression(model, P[t]- sum(β[k] * X[2*t,k] for k = K_d+1:K_d+K_s-1) + log(1 - θ *(β[2] + β[3] * X[2*t, end]) ) ) )
    end

    for l = 1:L
        push!(g, @NLexpression(model, sum(Z[t,l] * r[t] for t = 1:2*T)))
    end

    if estimation_method[2] == :log_constraint
        for t = 1:T
            @NLconstraint(model, 0 <= 1 - θ *(β[2] + β[3] * X[2*t, end]))
        end
    end
    @NLobjective(model, Min, sum( g[l] *Ω[l,k] * g[k] for l = 1:L, k = 1:L))
    optimize!(model)
    
    α_hat = value.(β)[1:K_d]
    γ_hat = value.(β)[K_d+1:end]
    θ_hat = value.(θ)

    if  sum(1 .- θ_hat .*(α_hat[2] .+ α_hat[3] .* X_s[:,end]) .<= 0) == 0
        return α_hat, γ_hat, θ_hat, termination_status_code(termination_status(model))
    else 
        error("The estimation result violates the model assumption ")
    end 
end

function GMM_estimation_MPEC(T, Q, P, Z, Z_s, Z_d, X, X_s, X_d, parameter, estimation_method::Tuple{Symbol, Symbol, Symbol}, start_value, tol_level)
    
    """ 
    Estimate the demand and supply parameter given a market simultaneously
    The 
    
    """
    @unpack γ_0, γ_1, γ_2, γ_3, θ, start_θ, start_γ = parameter

    L = size(Z, 2)
    K_s = size(X_s, 2)
    K_d = size(X_d, 2)

    # The wight function for the GMM estimation
    Ω = inv(Z' * Z)/T

    if tol_level == :tight
        tol = 1e-15
        acceptable_tol = 1e-12
    elseif tol_level == :loose
        tol = 1e-6
        acceptable_tol = 1e-5
    end

    if start_value == :true
        start_β = [α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3]
        start_θ = θ

    elseif start_value == :random
        start_β = [α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3] .+ rand(Uniform(-10, 10), 8)
        start_θ = θ + rand(Uniform(-10, 1))
    end

    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "tol", tol)
    set_optimizer_attribute(model, "max_iter", 1000)
    set_optimizer_attribute(model, "acceptable_tol", acceptable_tol)
    #set_silent(model)
    @variable(model, β[k = 1:K_d+K_s-1])

    if estimation_method[3] == :theta_constraint
        @variable(model, 0 <= θ <= 1)
    else
        @variable(model, θ)
    end

    MC = Any[];
    for t = 1:T
        push!(MC, @NLexpression(model, (1 - θ * (β[2] + β[3] * X[2*t, end]))* P[t]))
    end

    r = Any[];
    g = Any[];
    for t =1:T
        push!(r, @NLexpression(model, P[t] - sum(β[k] * X[2*t-1,k] for k = 1:K_d) ))
        push!(r, @NLexpression(model, log(MC[t]) - sum(β[k] * X[2*t,k] for k = K_d+1:K_d+K_s-1)))
    end

    for l = 1:L
        push!(g, @NLexpression(model, sum(Z[t,l] * r[t] for t = 1:2*T)))
    end

    if estimation_method[2] == :log_constraint
        for t = 1:T
            @NLconstraint(model, 0 <= 1 - θ *(β[2] + β[3] * X[2*t, end]))
        end
    end
    @NLobjective(model, Min, sum( g[l] *Ω[l,k] * g[k] for l = 1:L, k = 1:L))
    optimize!(model)
    
    α_hat = value.(β)[1:K_d]
    γ_hat = value.(β)[K_d+1:end]
    θ_hat = value.(θ)

    if  sum(1 .- θ_hat .*(α_hat[2] .+ α_hat[3] .* X_s[:,end]) .<= 0) == 0
        return α_hat, γ_hat, θ_hat, termination_status_code(termination_status(model))
    else 
        error("The estimation result violates the model assumption ")
    end 
end

function GMM_estimation_MPEC_linear(T, Q, P, Z, Z_s, Z_d, X, X_s, X_d, parameter, estimation_method::Tuple{Symbol, Symbol, Symbol}, start_value, tol_level)
    
    """ 
    Estimate the demand and supply parameter given a market simultaneously
    The 
    
    """
    @unpack γ_0, γ_1, γ_2, γ_3, θ, start_θ, start_γ = parameter

    L = size(Z, 2)
    K_s = size(X_s, 2)
    K_d = size(X_d, 2)

    # The wight function for the GMM estimation
    Ω = inv(Z' * Z)/T

    if tol_level == :tight
        tol = 1e-15
        acceptable_tol = 1e-12
    elseif tol_level == :loose
        tol = 1e-6
        acceptable_tol = 1e-5
    end

    if start_value == :true
        start_β = [α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3]
        start_θ = θ

    elseif start_value == :random
        start_β = [α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3] .+ rand(Uniform(-10, 10), 8)
        start_θ = θ + rand(Uniform(-10, 1))
    end

    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "tol", tol)
    set_optimizer_attribute(model, "max_iter", 1000)
    set_optimizer_attribute(model, "acceptable_tol", acceptable_tol)
    set_silent(model)
    @variable(model, β[k = 1:K_d+K_s-1])
    @variable(model, θ)

    MC = Any[];
    for t = 1:T
        push!(MC, @NLexpression(model,  P[t] + θ * (β[2] + β[3] * X[2*t, end])))
    end

    r = Any[];
    g = Any[];
    for t =1:T
        push!(r, @NLexpression(model, P[t] - sum(β[k] * X[2*t-1,k] for k = 1:K_d) ))
        push!(r, @NLexpression(model, MC[t] - sum(β[k] * X[2*t,k] for k = K_d+1:K_d+K_s-1)))
    end

    for l = 1:L
        push!(g, @NLexpression(model, sum(Z[t,l] * r[t] for t = 1:2*T)))
    end

    @NLobjective(model, Min, sum( g[l] *Ω[l,k] * g[k] for l = 1:L, k = 1:L))
    optimize!(model)
    
    α_hat = value.(β)[1:K_d]
    γ_hat = value.(β)[K_d+1:end]
    θ_hat = value.(θ)

    return α_hat, γ_hat, θ_hat, termination_status_code(termination_status(model))
end


function estimate_nonlinear_2SLS(parameter, data, estimation_method::Tuple{Symbol, Symbol, Symbol}, start_value, tol_level)
    """
    Given data, reshape the data and pass it to the function that implement the GMM estimation

    """

    @unpack T = parameter
    if estimation_method[1] == :mpec_linear
        Q  = data.Q
        w  = data.w
        r  = data.r
        p  = data.P
        y  = data.y
    else
        Q  = data.logQ
        w  = data.logw
        r  = data.logr
        p  = data.logP
        y  = data.logy
    end

    
    z  = data.z
    iv_w = data.iv_w
    iv_r = data.iv_r
    
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
        Z_st = vcat(1, z[t], w[t], r[t], y[t])
        Z_t = [Z_dt zeros(length(Z_dt));  zeros(length(Z_st)) Z_st]'

        X_dt = vcat(1, -Q[t], -z[t].*Q[t], y[t])
        X_st = vcat(1, Q[t], w[t], r[t], z[t])
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
        α_hat, γ_hat, θ_hat, status = GMM_estimation_separate(T, Q, P, Z, Z_s, Z_d, X, X_s, X_d, parameter, estimation_method, start_value, tol_level)
    elseif estimation_method[1] == :simultaneous
        α_hat, γ_hat, θ_hat, status = GMM_estimation_simultaneous(T, Q, P, Z, Z_s, Z_d, X, X_s, X_d, parameter, estimation_method , start_value, tol_level)
    elseif estimation_method[1] == :mpec
        α_hat, γ_hat, θ_hat, status = GMM_estimation_MPEC(T, Q, P, Z, Z_s, Z_d, X, X_s, X_d, parameter, estimation_method , start_value, tol_level)
    else
        α_hat, γ_hat, θ_hat, status = GMM_estimation_MPEC_linear(T, Q, P, Z, Z_s, Z_d, X, X_s, X_d, parameter, estimation_method , start_value, tol_level)
    end

    return α_hat, γ_hat, θ_hat, status
end

function iterate_esimation_nonlinear_2SLS(parameter, data, estimation_method::Tuple{Symbol, Symbol, Symbol}, start_value, tol_level)

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
        α_est_s, γ_est_s, θ_est_s, status_s = estimate_nonlinear_2SLS(parameter, data_s, estimation_method, start_value, tol_level)

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

function generate_contour_set_of_GMM(parameter, data, theta_range, gamma_range)
    """
    Draw the contour plot of the GMM value function with respect to the conduct parameter and the constant in the marginal cost function
    """

    @unpack α_0, α_1, α_2, α_3,γ_0 , γ_1 ,γ_2 ,γ_3, θ, σ ,T = parameter

    γ = [γ_0 , γ_1 ,γ_2 ,γ_3]
    α = [α_0, α_1, α_2, α_3]

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

    L = size(Z, 2)
    L_d = size(Z_d,2)
    L_s = size(Z_s,2)
    K_s = size(X_s, 2)
    K_d = size(X_d, 2)

    Ω = inv(Z_s' * Z_s)/T

    GMM_value = []


    for θ_hat = theta_range, γ_hat = gamma_range

        r = P .- γ_hat .-sum(γ[k] .* X_s[:,k] for k = 2:K_s-1) + log.(1 .- θ_hat .*(α_1 .+ α_2 .* X_s[:, end]) ) 

        g = sum(Z_s[t,:] .* r[t] for t = 1:T)

        gmm_value = sum( g[l] *Ω[l,k] * g[k] for l = 1:L_s, k = 1:L_s)

        push!(GMM_value, gmm_value)
    end

    return reshape(GMM_value, (length(gamma_range), length(theta_range)))
end

function generate_difference_GMM_value(parameter, data, estimation_result)
    """
    Draw the figure of the difference of the values of the GMM objective function under the true parameters and the estimated parameter  

    """
    
    @unpack α_0, α_1, α_2, α_3, γ_0, γ_1 ,γ_2 ,γ_3, θ, σ ,T = parameter

    γ = [γ_0 , γ_1 ,γ_2 ,γ_3]
    α = [α_0, α_1, α_2, α_3]

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

    L = size(Z, 2)
    L_d = size(Z_d,2)
    L_s = size(Z_s,2)
    K_s = size(X_s, 2)
    K_d = size(X_d, 2)

    Ω = inv(Z_s' * Z_s)/T
    θ_hat = estimation_result.θ
    γ_hat = estimation_result.γ_0


    
    if 1 <= sum(1 .- θ_hat .*(α_1 .+ α_2 .* X_s[:, end]) .<= 0)
        return missing, missing, missing
    else
        
        r_true = P .- γ_0 .-sum(γ[k] .* X_s[:,k] for k = 2:K_s-1) + log.(1 .- θ .*(α_1 .+ α_2 .* X_s[:, end]) ) 
        r_estimated = P .- γ_hat .-sum(γ[k] .* X_s[:,k] for k = 2:K_s-1) + log.(1 .- θ_hat .*(α_1 .+ α_2 .* X_s[:, end]) ) 

        g_true      = sum(Z_s[t,:] .* r_true[t] for t = 1:T)
        g_estimated = sum(Z_s[t,:] .* r_estimated[t] for t = 1:T)

        gmm_value_true = sum( g_true[l] *Ω[l,k] * g_true[k] for l = 1:L_s, k = 1:L_s)
        gmm_value_estimated = sum( g_estimated[l] *Ω[l,k] * g_estimated[k] for l = 1:L_s, k = 1:L_s)

        difference_theta = abs(θ_hat - θ)
        difference_gamma = abs(γ_hat - γ_0)
        difference_gmm_value = abs(gmm_value_true - gmm_value_estimated)

        return difference_theta, difference_gamma, difference_gmm_value
    end
end