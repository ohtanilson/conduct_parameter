using LinearAlgebra, Distributions
using Statistics, Random, MultivariateStats
using JuMP, Ipopt
using DelimitedFiles, JLD, CSV, DataFrames
using Plots, Combinatorics, Dates, StatsPlots
using Parameters: @unpack, @with_kw
using GLM


#---------------------------------------------------------------------------------------------------------

"""

Generate a simulation data under log demand and log marginal cost model

"""

market_parameters_log = @with_kw (
    α_0 = 10, # Demand parameter
    α_1 = 1,
    α_2 = 1,
    α_3 = 1,
    γ_0 = 1,  # Marginal cost parameter
    γ_1 = 1,
    γ_2 = 1,
    γ_3 = 1,
    θ = 0.5,  # Conduct paramter
    σ = 1,    # Standard deviation of the error term
    T = 50,   # Number of markets
    S = 1000, # Number of simulation
)

mutable struct market_data_log
    Q::Vector{Float64}
    W::Vector{Float64}
    R::Vector{Float64}
    Z::Vector{Float64}
    P::Vector{Float64}
    Y::Vector{Float64}
    MC::Vector{Float64}
    IV::Matrix{Float64}
end


function simulation_data_log(parameter)


    @unpack α_0, α_1, α_2, α_3,γ_0 , γ_1 ,γ_2 ,γ_3, θ,σ ,T, S = parameter
    
    Q = Float64[];
    P = Float64[];
    W = Float64[];
    R = Float64[];
    Z = Float64[];
    Y = Float64[];
    IV_W = Float64[];
    IV_R = Float64[];

    group_id = Int64[];
    for s = 1:S

        Random.seed!(s * 12345)
        for t = 1:T
            r_t = rand(Uniform(0,1))
            w_t = rand(Uniform(1,3))
            z_t = rand(Uniform(0,1))
            iv_w_t = w_t + randn()
            iv_r_t = r_t + randn()

            y_t = randn()
            ε_c = rand(Normal(0,σ))
            ε_d = rand(Normal(0,σ))

            log_Q_t = (α_0 + log(1 - θ * (α_1 + α_2 * z_t))  + α_3 * y_t - ( γ_0 + γ_2 * log(w_t) + γ_3 * log(r_t)) + ε_d- ε_c)/ (γ_1 + α_1 + α_2 * z_t) # Equilibrium total quantity

            log_p_t = α_0 - (α_1 + α_2 * z_t )* log_Q_t  + α_3 * y_t + ε_d # The demand function

            push!(group_id, s)
            push!(Q, log_Q_t)
            push!(P,log_p_t)
            push!(W,w_t)
            push!(R,r_t)
            push!(Z, z_t)
            push!(Y, y_t)
            push!(IV_W, iv_w_t)
            push!(IV_R, iv_r_t)

        end
    end

    # Save the data as DataFrame
    data = DataFrame(group_id_k = group_id, logQ = Q, logP = P, w = W, r = R, y = Y, z = Z, iv_w = IV_W, iv_r = IV_R)
    return data
end



function simulation_data_log_without_demand_shifter(parameter)


    @unpack α_0, α_1, α_2,γ_0 , γ_1 ,γ_2 ,γ_3, θ,σ ,T, S = parameter
    
    Q = Float64[];
    P = Float64[];
    W = Float64[];
    R = Float64[];
    Z = Float64[];
    IV_W = Float64[];
    IV_R = Float64[];

    group_id = Int64[];
    for s = 1:S

        Random.seed!(s * 12345)
        for t = 1:T
            r_t = rand(Uniform(0,1))
            z_t = rand(Uniform(0,1))
            w_t = rand(Uniform(1,3))
            iv_w_t = w_t + randn()
            iv_r_t = r_t + randn()

            ε_c = rand(Normal(0,σ))
            ε_d = rand(Normal(0,σ))

            log_Q_t = (α_0 + log(1 - θ * (α_1 + α_2 * z_t))  - ( γ_0 + γ_2 * log(w_t) + γ_3 * log(r_t)) + ε_d- ε_c)/ (γ_1 + α_1 + α_2 * z_t) # Equilibrium total quantity

            log_p_t = α_0 - (α_1 + α_2 * z_t )* log_Q_t + ε_d # The demand function

            push!(group_id, s)
            push!(Q, log_Q_t)
            push!(P,log_p_t)
            push!(W,w_t)
            push!(R,r_t)
            push!(Z, z_t)
            push!(IV_W, iv_w_t)
            push!(IV_R, iv_r_t)

        end
    end

    # Save the data as DataFrame
    data = DataFrame(group_id_k = group_id, logQ = Q, logP = P, w = W, r = R, z = Z, iv_w = IV_W, iv_r = IV_R)
    return data
end



#---------------------------------------------------------------------------------------------------------------------

# Plotting logP and logQ

for theta = [0.1:0.1:0.9;]

    parameter_plot = market_parameters_log(θ = theta, α_2 = 0.1);

    data_plot = simulation_data_log(parameter_plot);

    # Plot the scatter plot of P and Q

    ols_plot = lm(@formula(logP ~ 1 + logQ), data_plot[1:500,:])
    plot_logP_logQ = plot(scatter(data_plot.logP[1:500], data_plot.logQ[1:500], ms=2, ma=0.2), ylabel = "logP", xlabel = "logQ", aspect_ratio=:equal, title = " θ =  $theta")
    plot!(data_plot.logQ[1:500], predict(ols_plot), label = "OLS")

    savefig(plot_logP_logQ, "../conduct_parameter/Yuri/plot_theta_$theta.pdf")

end







# Generate the simulation data and save the data as CSV file
# Simulation data wit demand shifter y  
for t in [50, 100, 200, 1000],  sigma in [0.001, 0.5, 1, 2] 
        
    parameter_plot = market_parameters_log(T = t, σ = sigma, α_2 = 0.1);
    data_plot = simulation_data_log(parameter_plot);


    ols_plot = lm(@formula(logP ~ 1 + logQ), data_plot[1:(t*10),:])
    plot_logP_logQ = plot(scatter(data_plot.logP[1:t*10], data_plot.logQ[1:t*10], ms=2, ma=0.2), ylabel = "logP", xlabel = "logQ", aspect_ratio=:equal, title = " t =  $t, σ = $sigma")
    plot!(data_plot.logQ[1:t*10], predict(ols_plot), label = "OLS")

    filename = "../conduct_parameter/Yuri/plot_"*string(t)*"_and_"*string(sigma)*".pdf"

    savefig(plot_logP_logQ, filename)
    
end





#---------------------------------------------------------------------------------------------------------------------





# Generate the simulation data and save the data as CSV file
# Simulation data wit demand shifter y  
for t in [50, 100, 200, 1000],  sigma in [0.001, 0.5, 1, 2] 
        
    parameter = market_parameters_log(T = t, σ = sigma);
    data = simulation_data_log(parameter);


    if sigma == 1 || sigma == 2
        sigma = Int64(sigma)
    end

    @unpack α_1, α_2, θ = parameter;

    # Save the simulation data if the demand function is downward sloping
    if sum(0 .<= - (α_1 .+ α_2 * data.z)) == 0 && sum(1 .- θ * (α_1 .+ α_2 * data.z) .<= 0) == 0 

        file_name = "../conduct_parameter/output/data_loglinear_loglinear_n_"*string(t)*"_sigma_"*string(sigma)*".csv"

        CSV.write(file_name, data)
    else
        error("The demand function is not downward sloping")
    end
end


# Generate the simulation data and save the data as CSV file
# Simulation data without demand shifter y    
for t in [50, 100, 200, 1000],  sigma in [0.001, 0.5, 1, 2] 
        
    parameter = market_parameters_log(T = t, σ = sigma, α_3 = 0);
    data = simulation_data_log(parameter);


    if sigma == 1 || sigma == 2
        sigma = Int64(sigma)
    end

    @unpack α_1, α_2, θ = parameter;

    # Save the simulation data if the demand function is downward sloping
    if sum(0 .<= - (α_1 .+ α_2 * data.z)) == 0 && sum(1 .- θ * (α_1 .+ α_2 * data.z) .<= 0) == 0 

        file_name = "../conduct_parameter/output/data_loglinear_loglinear_n_"*string(t)*"_sigma_"*string(sigma)*"_without_demand_shifter_y"*".csv"

        CSV.write(file_name, data)
    else
        error("The demand function is not downward sloping")
    end
end