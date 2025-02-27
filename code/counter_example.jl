using Plots
using LinearAlgebra  # 線形代数の計算用
using Statistics     # 統計計算用

# Define the functions
P(Q, X) = X * exp(-Q) + 1
∂P∂Q(Q, X) = -Q * X * exp(-Q)  # This is ∂P/∂Q × Q

# Test points
Q1, X1 = 1.0, 2.0
Q2, X2 = 2.0, 2exp(1)  # ≈ 5.437

# Calculate P values
P1 = P(Q1, X1)
P2 = P(Q2, X2)

# Calculate partial derivatives
dP1 = ∂P∂Q(Q1, X1)
dP2 = ∂P∂Q(Q2, X2)

# Print results
println("Comparing two points:")
println("Point 1: (Q=$Q1, X=$X1)")
println("Point 2: (Q=$Q2, X=$X2)")
println("\nP values:")
println("P(Q1,X1) = $P1")
println("P(Q2,X2) = $P2")
println("Difference in P values: $(abs(P1-P2))")
println("\nPartial derivatives:")
println("∂P/∂Q × Q at point 1: $dP1")
println("∂P/∂Q × Q at point 2: $dP2")
println("Difference in derivatives: $(abs(dP1-dP2))")

## Visualization

Q_range = 0:0.1:5
X_values = [2.0, 2exp(1)]

# Plot P(Q,X) for different X values
plot_counter_example = plot(Q_range, Q -> P(Q, X1), label="G(Q_1,X_1;0)", 
     xlabel="Q", ylabel="Value", 
     title="P(Q,X) = X exp(-Q) + 1", color = :blue)
plot!(Q_range, Q -> P(Q, X2), label="G(Q_2,X_2;0)", color = :red)
hline!([P1], linestyle=:dash, color=:black, label="")
hline!([0], linestyle=:dash, color=:gray, label="", alpha=0.5)


# Plot ∂P/∂Q × Q for different X values
plot!(Q_range, Q -> ∂P∂Q(Q, X1)+ P(Q, X1), label="G(Q_1,X_1;1)", linestyle=:dash, color = :blue)
plot!(Q_range, Q -> ∂P∂Q(Q, X2)+ P(Q, X2), label="G(Q_2,X_2;1)", linestyle=:dash, color = :red)

# Add points of interest
scatter!([Q1, Q2], [P1, P2], label="", color=:black)




# Save the plot
savefig(plot_counter_example, "../conduct_parameter/figuretable/plot_counter_example.pdf")




########################################################

## Specific function in Theorem 

########################################################

θ = 1/2
θ_alternative = 1/3

# Define the functions
P(Q, X) = Q^(-1/θ) * X 
∂P∂Q(Q, X, θ) = (-1/θ) * Q^(-1/θ-1) * X 

G(Q, X, θ, θ_alternative) = P(Q, X) + θ_alternative * ∂P∂Q(Q, X, θ) * Q


Q_range = 0.5:0.01:1
X_values = [0.1, 0.11]

# Plot P(Q,X) for different X values
plot_counter_example = plot(Q_range, Q -> G(Q, X1, θ, θ), label="G(Q_1,X_1;θ)", 
     xlabel="Q", ylabel="Value", 
     title="P(Q,X) = X exp(-Q) + 1", color = :blue)
plot!(Q_range, Q -> G(Q, X2, θ, θ), label="G(Q_2,X_2;θ)", color = :red)
hline!([0], linestyle=:dash, color=:gray, label="", alpha=0.5)


# Plot ∂P/∂Q × Q for different X values
plot!(Q_range, Q -> G(Q, X1, θ, θ_alternative), label="G(Q_1,X_1;θ)", linestyle=:dash, color = :blue)
plot!(Q_range, Q -> G(Q, X2, θ, θ_alternative), label="G(Q_2,X_2;θ)", linestyle=:dash, color = :red)







