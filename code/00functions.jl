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

function GMM_estimation_separate(T, Q, P, Z, Z_s, Z_d, X, X_s, X_d, α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3, θ_0, start_θ, start_γ, estimation_method::Tuple{Symbol, Symbol, Symbol}, starting_value, tol_level)

    """
    Estimate the demand and supply parameter given a market
    The demand parameters are destimated by IV.
    The supply parameters and the conduct parameter are estimated by the GMM.
    """
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


    start_γ = zeros(4)
    start_θ = 0
    if starting_value == :true_value
        start_θ = θ_0
        start_γ = [γ_0, γ_1, γ_2, γ_3]

    elseif starting_value == :random
        start_γ = [γ_0, γ_1, γ_2, γ_3] .+ rand(Uniform(-20, 20), 4)
        start_θ = 10
        while sum(1 .- start_θ .*(α_hat[2] .+ α_hat[3] .* X_s[:,end]) .<= 0) != 0
            start_θ = θ_0 + rand(Uniform(-10, 1))
        end
    end

    # Pick up the index of market under which the inside of the log has a negative value
    for t = 1:T
        if 1 .- θ_0 .*(α_hat[2] .+ α_hat[3] .* X_s[t,end]) <= 0
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
        return α_hat, γ_hat, θ_hat, termination_status_code(JuMP.termination_status(model))
        #if  sum(1 .- θ_hat .*(α_hat[2] .+ α_hat[3] .* X_s[:,end]) .<= 0) == 0

        #   return α_hat, γ_hat, θ_hat, termination_status_code(JuMP.termination_status(model))
        #else
        #    error("The estimation result violates the model assumption ")
        #end
    end

end


function GMM_estimation_simultaneous(T, Q, P, Z, Z_s, Z_d, X, X_s, X_d, α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3, θ_0, start_θ, start_γ, estimation_method::Tuple{Symbol, Symbol, Symbol}, starting_value, tol_level)
    
    """ 
    Estimate the demand and supply parameter given a market simultaneously
    The 
    
    """

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

    start_β = zeros(8)
    start_θ = 0

    if starting_value == :true_value
        start_β = [α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3]
        start_θ = θ_0

    elseif starting_value == :random
        start_β = [α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3] .+ rand(Uniform(-10, 10), 8)
        while sum(1 .- start_θ .*(α_1[2] .+ α_2 .* X[:,end]) .<= 0) != 0
            start_θ = θ_0 + rand(Uniform(-10, 1))
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

    return α_hat, γ_hat, θ_hat, termination_status_code(termination_status(model))
    #if  sum(1 .- θ_hat .*(α_hat[2] .+ α_hat[3] .* X_s[:,end]) .<= 0) == 0
    #    return α_hat, γ_hat, θ_hat, termination_status_code(termination_status(model))
    #else 
    #    error("The estimation result violates the model assumption ")
    #end 
end

function GMM_estimation_MPEC(T, Q, P, Z, Z_s, Z_d, X, X_s, X_d, α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3, θ_0, start_θ, start_γ,estimation_method::Tuple{Symbol, Symbol, Symbol}, starting_value, tol_level)
    
    """ 
    Estimate the demand and supply parameter given a market simultaneously

    """

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

    start_β = zeros(8)
    start_θ = 0

    if starting_value == :true_value
        start_β = [α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3]
        start_θ = θ_0
    elseif starting_value == :random
        start_β = [α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3] .+ rand(Uniform(-10, 10), 8)
        start_θ = θ_0 + rand(Uniform(-10, 1))
    end


    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "tol", tol)
    set_optimizer_attribute(model, "max_iter", 1000)
    set_optimizer_attribute(model, "acceptable_tol", acceptable_tol)
    set_silent(model)
    @variable(model, β[k = 1:K_d+K_s-1], start = start_β[k])

    #@constraint(model, c1, β[1] >=0) # constant term should be positive
    #@constraint(model, c2, β[K_d+1] >=0) # constant term should be positive
    #@constraint(model, c3, β[2] >=0) # demand curve should be downward
    #@constraint(model, c4, β[3] >=0) # demand curve should be downward

    if estimation_method[3] == :theta_constraint
        @variable(model, 0 <= θ <= 1, start = start_θ)
    else
        @variable(model, θ, start = start_θ)
    end

    @variable(model, 0 <= MC[t = 1:T])
    for t = 1:T
        @NLconstraint(model, exp(P[t]) == MC[t] + θ * (β[2] + β[3] * X[2*t, end])* exp(P[t]))
    end

    # MC = Any[];
    # for t = 1:T
    #    push!(MC, @NLexpression(model, (1 - θ * (β[2] + β[3] * X[2*t, end]))* exp(P[t])))
    # end
    # if estimation_method[2] == :log_constraint
    #    for t = 1:T
    #        @NLconstraint(model, 0 <= 1 - θ *(β[2] + β[3] * X[2*t, end]))
    #    end
    # end

    r = Any[];
    g = Any[];
    for t =1:T
        push!(r, @NLexpression(model, P[t] - sum(β[k] * X[2*t-1,k] for k = 1:K_d) ))
        push!(r, @NLexpression(model, log(MC[t]) - sum(β[k] * X[2*t,k] for k = K_d+1:K_d+K_s-1)))
    end

    for l = 1:L
        push!(g, @NLexpression(model, sum(Z[t,l] * r[t] for t = 1:2*T)))
    end
    @NLobjective(model, Min, sum( g[l] * Ω[l,k] * g[k] for l = 1:L, k = 1:L))
    optimize!(model)
    
    α_hat = value.(β)[1:K_d]
    γ_hat = value.(β)[K_d+1:end]
    θ_hat = value.(θ)

    return α_hat, γ_hat, θ_hat, termination_status_code(termination_status(model))
    #if  sum(1 .- θ_hat .*(α_hat[2] .+ α_hat[3] .* X_s[:,end]) .<= 0) == 0
    #    return α_hat, γ_hat, θ_hat, termination_status_code(termination_status(model))
    #else 
    #    error("The estimation result violates the model assumption ")
    #end 
end

function GMM_estimation_Optim(T, Q, P, Z, Z_s, Z_d, X, X_s, X_d, α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3, θ_0, start_θ, start_γ,estimation_method::Tuple{Symbol, Symbol, Symbol}, starting_value, tol_level)
    
    """ 
    Estimate the demand and supply parameter given a market simultaneously
    """
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

    #if starting_value == :true_value
        global start_β = [α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3]
        global start_θ = θ_0

    #elseif starting_value == :random
        global start_β = [α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3] .+ rand(Uniform(-10, 10), 8)
        global start_θ = θ_0 + rand(Uniform(-10, 1))
    #end

    # model = Model(Ipopt.Optimizer)
    # set_optimizer_attribute(model, "tol", tol)
    # set_optimizer_attribute(model, "max_iter", 1000)
    # set_optimizer_attribute(model, "acceptable_tol", acceptable_tol)
    # set_silent(model)
    #@variable(model, β[k = 1:K_d+K_s-1])
    # if estimation_method[3] == :theta_constraint
    #     θ = zeros(1)

    # else
    
    # end
    
    
    f = function(target_param)
        β = target_param[1:K_d+K_s-1]
        if estimation_method[3] == :theta_constraint
            # logit: 0 <= x/(1+x) <= 1
            θ = target_param[K_d+K_s]/(1+target_param[K_d+K_s])
        else
            θ = target_param[K_d+K_s]
        end
        β[1] = exp(β[1]) # constant term should be positive
        β[K_d+1] = exp(β[K_d+1]) # constant term should be positive
        β[2] = exp(β[2]) # demand curve should be downward
        β[3] = exp(β[3]) # demand curve should be downward
    
        
        MC = zeros(T);
        for t = 1:T
            MC[t] = (1.0 - θ * (β[2] + β[3] * X[2*t, end]))* exp(P[t])
        end
        r = zeros(2*T);
        g = zeros(T);
        for t = 1:T
            xb = 0
            for k = 1:K_d
                xb += β[k] * X[2*t-1,k]    
            end
            r[2*t-1] = P[t] - xb
            for k = K_d+1:K_d+K_s-1
                xb += β[k] * X[2*t,k]  
            end
            r[2*t] = log(MC[t]) - xb
        end
        for l = 1:L
            Zr = 0
            for t = 1:2*T
                Zr += Z[t,l] * r[t]
            end
            g[l] = Zr
        end
        objQ = 0
        for l = 1:L, k = 1:L
            objQ += g[l] *Ω[l,k] * g[k]
        end
        return objQ
    end

    initial_x = vcat(start_β,start_θ)
    results = Optim.optimize(f, initial_x)
    status = Optim.converged(results)
    α_hat = results.minimizer[1:K_d]
    α_hat[1] = exp(α_hat[1]) # constant term should be positive
    α_hat[2] = exp(α_hat[2]) # demand curve should be downward
    α_hat[3] = exp(α_hat[3]) # demand curve should be downward
    γ_hat = results.minimizer[K_d+1:K_d+K_s-1]
    γ_hat[1] = exp(γ_hat) # constant term should be positive


    if estimation_method[3] == :theta_constraint
        # logit: 0 <= x/(1+x) <= 1
        θ_hat = results.minimizer[K_d+K_s]/(1+results.minimizer[K_d+K_s])
    else
        θ_hat = results.minimizer[K_d+K_s]
    end

    return α_hat, γ_hat, θ_hat, termination_status_code(termination_status(model))
    #if  sum(1 .- θ_hat .*(α_hat[2] .+ α_hat[3] .* X_s[:,end]) .<= 0) == 0
    #    return α_hat, γ_hat, θ_hat, status#termination_status_code(termination_status(model))
    #else 
    #    error("The estimation result violates the model assumption ")
    #end 
end

function GMM_estimation_MPEC_linear(T, Q, P, Z, Z_s, Z_d, X, X_s, X_d, α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3, θ_0, start_θ, start_γ,  estimation_method::Tuple{Symbol, Symbol, Symbol}, starting_value, tol_level)
    
    """ 
    Estimate the demand and supply parameter given a market simultaneously
    """

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

    if starting_value == :true_value
        start_β = [α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3]
        start_θ = θ

    elseif starting_value == :random
        start_β = [α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3] .+ rand(Uniform(-10, 10), 8)
        start_θ = θ + rand(Uniform(-10, 1))
    end

    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "tol", tol)
    set_optimizer_attribute(model, "max_iter", 1000)
    set_optimizer_attribute(model, "acceptable_tol", acceptable_tol)
    set_silent(model)
    @variable(model, β[k = 1:K_d+K_s-1])
    @variable(model, 0 <= θ <= 1)

    MC = Any[];
    for t = 1:T
        push!(MC, @NLexpression(model,  P[t] - θ * (β[2] + β[3] * X[2*t, end]) * X[2*t,K_d + 2]   ))
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



@everywhere function estimate_nonlinear_2SLS(simulation_setting::SIMULATION_SETTING)
    """
    Given data, reshape the data and pass it to the function that implement the GMM estimation

    """

    α_0 = simulation_setting.α_0
    α_1 = simulation_setting.α_1
    α_2 = simulation_setting.α_2
    α_3 = simulation_setting.α_3
    γ_0 = simulation_setting.γ_0
    γ_1 = simulation_setting.γ_1
    γ_2 = simulation_setting.γ_2
    γ_3 = simulation_setting.γ_3
    θ_0 = simulation_setting.θ_0
    start_θ = simulation_setting.start_θ
    start_γ = simulation_setting.start_γ
    T = simulation_setting.T
    data = simulation_setting.data
    estimation_method = simulation_setting.estimation_method
    starting_value = simulation_setting.starting_value
    tol_level = simulation_setting.tol_level
    simulation_index  = simulation_setting.simulation_index

    data_s = data[(simulation_index-1)*T+1:simulation_index*T,:]

    if estimation_method[1] == :mpec_linear
        Q  = data_s.Q
        w  = data_s.w
        r  = data_s.r
        p  = data_s.P
        y  = data_s.y
    else
        Q  = data_s.logQ
        w  = data_s.logw
        r  = data_s.logr
        p  = data_s.logP
        y  = data_s.logy
    end

    
    z  = data_s.z
    iv_w = data_s.iv_w
    iv_r = data_s.iv_r
    
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
        α_hat, γ_hat, θ_hat, status = GMM_estimation_separate(T, Q, P, Z, Z_s, Z_d, X, X_s, X_d, α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3, θ_0, start_θ, start_γ , estimation_method, starting_value, tol_level)
    elseif estimation_method[1] == :simultaneous
        α_hat, γ_hat, θ_hat, status = GMM_estimation_simultaneous(T, Q, P, Z, Z_s, Z_d, X, X_s, X_d, α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3, θ_0, start_θ, start_γ, estimation_method , starting_value, tol_level)
    elseif estimation_method[1] == :mpec
        α_hat, γ_hat, θ_hat, status = GMM_estimation_MPEC(T, Q, P, Z, Z_s, Z_d, X, X_s, X_d,  α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3, θ_0, start_θ, start_γ ,estimation_method , starting_value, tol_level)
    elseif estimation_method[1] == :mpec_linear
        α_hat, γ_hat, θ_hat, status = GMM_estimation_MPEC_linear(T, Q, P, Z, Z_s, Z_d, X, X_s, X_d, α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3, θ_0, start_θ, start_γ, estimation_method , starting_value, tol_level)
    elseif estimation_method[1] == :optim_nelder_mead
        α_hat, γ_hat, θ_hat, status = GMM_estimation_Optim(T, Q, P, Z, Z_s, Z_d, X, X_s, X_d, α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3, θ_0, start_θ, start_γ, estimation_method , starting_value, tol_level)
    end

    return α_hat, γ_hat, θ_hat, status
end

@everywhere function iterate_esimation_nonlinear_2SLS(parameter, data, estimation_method::Tuple{Symbol, Symbol, Symbol}, starting_value, tol_level)

    """
    Given the simulation data, run the estimation in each simulation index s = 1,..., 1000, and store the simulation results as a DataFrame file.
    """

    @unpack α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3, θ_0, start_θ, start_γ, T, S, σ = parameter

    α_est  = Vector{Float64}[]
    γ_est  = Vector{Union{Missing, Float64}}[]
    θ_est  = Union{Missing, Float64}[]
    status = Union{Missing, Int64}[]

    simulation_setting = [SIMULATION_SETTING(α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3, θ_0, start_θ, start_γ, T, σ, data, estimation_method, starting_value, tol_level,simulation_index) for simulation_index = 1:S]

    @elapsed @time simulation_mpec_result = pmap(estimate_nonlinear_2SLS, simulation_setting);

    for s = 1:S
        push!(α_est, simulation_mpec_result[s][1])
        push!(γ_est, simulation_mpec_result[s][2])
        push!(θ_est, simulation_mpec_result[s][3])
        push!(status,simulation_mpec_result[s][4])
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

    @unpack α_0, α_1, α_2, α_3,γ_0 , γ_1 ,γ_2 ,γ_3, θ_0, σ ,T = parameter

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
    
    @unpack α_0, α_1, α_2, α_3, γ_0, γ_1 ,γ_2 ,γ_3, θ_0, σ ,T = parameter

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