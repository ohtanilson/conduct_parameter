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

function GMM_estimation_linear_separate(T, P, X_s, X_d, Z_d, Z_s, Ω, γ_0, γ_1, γ_2, γ_3, θ_0, starting_value, tol_level)


    """
    Estimate the demand and supply parameter given a market
    The demand parameters are destimated by IV.
    The supply parameters and the conduct parameter are estimated by the GMM.
    """

    # first stage
    QZ_hat = Z_d * inv(Z_d' * Z_d) * Z_d' * (Z_d[:,2] .* X_d[:,2])
    Q_hat = Z_d * inv(Z_d' * Z_d) * Z_d' *  X_d[:,2]
    # second stage
    X_dd = hcat(ones(T), -Q_hat, -QZ_hat, X_d[:,end])

    α_hat = inv(X_dd' * X_dd) * (X_dd' * P)

    L_s = size(Z_s,2)
    K_s = size(X_s,2) - 1

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
        start_γ = [γ_0, γ_1, γ_2, γ_3] .+ rand(Uniform(-10, 10), 4)
        start_θ = θ_0 + rand(Uniform(-10, 10))

    end
    
    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "tol", tol)
    set_optimizer_attribute(model, "max_iter", 1000)
    set_optimizer_attribute(model, "acceptable_tol", acceptable_tol)
    set_silent(model)
    @variable(model, γ[k = 1:K_s], start = start_γ[k])
    @variable(model, 0 <= θ <= 1, start = start_θ)

    r = Any[];
    for t =1:T
        push!(r, JuMP.@NLexpression(model, P[t]- sum(γ[k] * X_s[t,k] for k = 1:K_s) - θ *(α_hat[2] + α_hat[3] * X_s[t, end]) * X_s[t, 2] )      )
    end

    g = Any[];
    for l = 1:L_s
        push!(g, JuMP.@NLexpression(model, sum(Z_s[t,l] * r[t] for t = 1:T)))
    end

    JuMP.@NLobjective(model, Min, sum(g[l]*Ω[l,k] * g[k] for l = 1:L_s, k = 1:L_s))
    JuMP.optimize!(model)

    γ_hat = JuMP.value.(γ)
    θ_hat = JuMP.value.(θ)

    ε_d = P .- sum(α_hat[k] * X_d[:, k] for k = 1:K_d)                                                    # T × 1 vector
    ε_s = P .- sum(γ_hat[k] * X_s[:, k] for k = 1:K_s) .- θ_hat * (α_hat[2] .+ α_hat[3] .* X_s[:,end]) .* X_s[:, 2]


    variance_demand = (X_d' * X_d)^(-1) * (sum(ε_d[t].^2 * X_d[t,:] * X_d[t,:]' for t = 1:T)/(T - K_d)) * (X_d' * X_d)^(-1)
    variacne_supply = (X_s' * X_s)^(-1) * (sum(ε_s[t].^2 * X_s[t,1:K_s-1] * X_s[t,1:K_s-1]' for t = 1:T))/(T - K_s) * (X_s' * X_s)^(-1)

    return α_hat, γ_hat, θ_hat, variance_demand, variacne_supply, termination_status_code(JuMP.termination_status(model))
end


function GMM_estimation_linear_simultaneous(T, P, Z, X, X_s, X_d, Ω, α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3, θ_0, starting_value, tol_level)
    
    """ 
    Estimate the demand and supply parameter given a market simultaneously
    """

    L = size(Z, 2)
    K_s = size(X_s, 2)
    K_d = size(X_d, 2)


    if tol_level == :tight
        tol = 1e-15
        acceptable_tol = 1e-12
    elseif tol_level == :loose
        tol = 1e-6
        acceptable_tol = 1e-5
    end

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
    @variable(model, β[k = 1:K_d+K_s-1])
    @variable(model, 0 <= θ <= 1)


    r = Any[];
    for t =1:T
        push!(r, @NLexpression(model, P[t] - sum(β[k] * X[2*t-1,k] for k = 1:K_d) ))
        push!(r, @NLexpression(model, P[t] - θ * (β[2] + β[3] * X[2*t, end]) * X[2*t, K_d+2] - sum(β[k] * X[2*t,k] for k = K_d+1:K_d+K_s-1)))
    end

    g = Any[];
    for l = 1:L
        push!(g, @NLexpression(model, sum(Z[t,l] * r[t] for t = 1:2*T)))
    end

    @NLobjective(model, Min, sum( g[l] *Ω[l,k] * g[k] for l = 1:L, k = 1:L))
    optimize!(model)
    
    α_hat = value.(β)[1:K_d]
    γ_hat = value.(β)[K_d+1:end]
    θ_hat = value.(θ)

    
    ε_d = P .- sum(α_hat[k] * X_d[:, k] for k = 1:K_d)                                                    # T × 1 vector
    ε_s = P .- sum(γ_hat[k] * X_s[:, k] for k = 1:K_s) .- θ_hat * (α_hat[2] .+ α_hat[3] .* X_s[:,end]) .* X_s[:, 2]


    variance_demand = (X_d' * X_d)^(-1) * (sum(ε_d[t].^2 * X_d[t,:] * X_d[t,:]' for t = 1:T)/(T - K_d)) * (X_d' * X_d)^(-1)
    variacne_supply = (X_s' * X_s)^(-1) * (sum(ε_s[t].^2 * X_s[t,1:K_s-1] * X_s[t,1:K_s-1]' for t = 1:T))/(T - K_s) * (X_s' * X_s)^(-1)

    return α_hat, γ_hat, θ_hat, variance_demand, variacne_supply, termination_status_code(JuMP.termination_status(model))
end




function compute_optimal_instruments(T, α_hat, γ_hat, θ_hat, X_d, X_s, P, Q_fitted_on_Z)

    K_d = size(X_d, 2)
    K_s = size(X_s, 2) - 1


    ε_d = P .- sum(α_hat[k] * X_d[:, k] for k = 1:K_d)                                                    # T × 1 vector
    ε_s = P .- sum(γ_hat[k] * X_s[:, k] for k = 1:K_s) .- θ_hat * (α_hat[2] .+ α_hat[3] .* X_s[:,end]) .* X_s[:, 2]

    ε = hcat(ε_d, ε_s)  # T × 2 mateix 


    Ω_estimate = inv(ε' * ε)/T


    Q_bar = Q_fitted_on_Z
    Q_bar_Z = X_s[:,end] .* Q_bar
    Q_bar_α_Z = (α_hat[2] .+ α_hat[3] .* X_s[:,end]) .* Q_bar

    Z_optimal = []

    for t = 1:T

        #Ω_t = inv(ε[t,:] * ε[t,:]') + diagm(ones(2) * 10e-8)

        D_z = [-1 Q_bar[t] Q_bar_Z[t] (-X_d[t,4]) 0 0 0 0 0; 0 (- θ_hat * Q_bar[t]) (θ_hat*Q_bar_Z[t]) 0 (-1)  (-Q_bar[t]) (-X_s[t,3]) (- X_s[t,4]) (-Q_bar_α_Z[t])]

        Z_t_optimal = Ω_estimate * D_z

        push!(Z_optimal, Z_t_optimal)
    end

    Z_optimal = reduce(vcat, Z_optimal)

    return Z_optimal
end


##################################################################################################################################################





##################################################################################################################################################


function estimation_linear_GMM_optimal_instrument(simulation_setting::SIMULATION_SETTING)
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


    Q  = data_s.Q
    w  = data_s.w
    r  = data_s.r
    p  = data_s.P
    y  = data_s.y

    z  = data_s.z
    iv_w = data_s.iv_w
    iv_r = data_s.iv_r

    q_fitted_on_Z = data_s.fitted_values_of_quantity_on_z
    
    iv = hcat(iv_w, iv_r)

    X_d = []
    X_s = []
    X   = []
    Z   = []
    Z_d = []
    Z_s = []
    P   = []
    Q_fitted_on_Z = []

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
        push!(Q_fitted_on_Z, q_fitted_on_Z[t])
    end

    Z   = reduce(vcat,(Z))
    X   = reduce(vcat,(X))
    X_d = reduce(vcat,(X_d))
    X_s = reduce(vcat,(X_s))
    Z_d = reduce(vcat,(Z_d))
    Z_s = reduce(vcat,(Z_s))


    if estimation_method == :linear_optimal_separate

        # The wight function for the GMM estimation
        Ω_initial = inv(Z_s' * Z_s)/T

        α_hat, γ_hat, θ_hat, variance_demand, variacne_supply,status = GMM_estimation_linear_separate(T, P, X_s, X_d, Z_d, Z_s, Ω_initial, γ_0, γ_1, γ_2, γ_3, θ_0, starting_value, tol_level)

        Z_optimal = compute_optimal_instruments(T, α_hat, γ_hat, θ_hat, X_d, X_s, P, Q_fitted_on_Z)
        Ω_optimal = inv(Z_optimal' * Z_optimal)/T
    
        α_optimal, γ_optimal, θ_optimal, variance_demand, variacne_supply, status = GMM_estimation_linear_separate(T, P, X_s, X_d, Z_d, Z_s, I_weight, γ_0, γ_1, γ_2, γ_3, θ_0, starting_value, tol_level)

        #return α_hat, γ_hat, θ_hat, variance_demand, variacne_supply, status
        return  α_optimal, γ_optimal, θ_optimal, variance_demand, variacne_supply, status
    
    elseif estimation_method == :linear_optimal_simultaneous

        # The wight function for the GMM estimation
        Ω_initial = inv(Z' * Z)/T

        α_hat, γ_hat, θ_hat, variance_demand, variacne_supply, status = GMM_estimation_linear_simultaneous(T, P, Z, X, X_s, X_d, Ω_initial, α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3, θ_0, starting_value, tol_level)

        Z_optimal = compute_optimal_instruments(T, α_hat, γ_hat, θ_hat, X_d, X_s, P, Q_fitted_on_Z)
        Ω_optimal = inv(Z_optimal' * Z_optimal)/T
    
        α_optimal, γ_optimal, θ_optimal, variance_demand, variacne_supply, status = GMM_estimation_linear_simultaneous(T, P, Z_optimal, X, X_s, X_d, Ω_optimal, α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3, θ_0, starting_value, tol_level)

        #return α_hat, γ_hat, θ_hat, status
        return  α_optimal, γ_optimal, θ_optimal, variance_demand, variacne_supply, status
    end
end

function simulation_GMM_optimal_instrument(parameter, data, estimation_method::Symbol, starting_value, tol_level)

    """
    Given the simulation data, run the estimation in each simulation index s = 1,..., 1000, and store the simulation results as a DataFrame file.
    """

    @unpack α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3, θ_0, start_θ, start_γ, T, S, σ = parameter

    α_est  = Vector{Float64}[]
    γ_est  = Vector{Union{Missing, Float64}}[]
    θ_est  = Union{Missing, Float64}[]
    variance_demand_est = Union{Missing, Float64}[]
    variacne_supply_est = Union{Missing, Float64}[]
    status = Union{Missing, Int64}[]

    simulation_setting = [SIMULATION_SETTING(α_0, α_1, α_2, α_3, γ_0, γ_1, γ_2, γ_3, θ_0, start_θ, start_γ, T, σ, data, estimation_method, starting_value, tol_level,simulation_index) for simulation_index = 1:S]

    @elapsed @time simulation_mpec_result = pmap(estimation_linear_GMM_optimal_instrument, simulation_setting);

    for s = 1:S
        push!(α_est, simulation_mpec_result[s][1])
        push!(γ_est, simulation_mpec_result[s][2])
        push!(θ_est, simulation_mpec_result[s][3])
        push!(variance_demand_est, simulation_mpec_result[s][4])
        push!(variacne_supply_est, simulation_mpec_result[s][5])
        push!(status,simulation_mpec_result[s][6])
    end 

    α_est = reduce(vcat, α_est')
    γ_est = reduce(vcat, γ_est')
    θ_est = reduce(vcat, θ_est')
    variance_demand_est = reduce(vcat, variance_demand_est)
    variacne_supply_est = reduce(vcat, variacne_supply_est)
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
    variance_demand = variance_demand_est,
    variacne_supply = variacne_supply_est,
    status = status,
    status_indicator = status_indicator)

    return estimation_result
end