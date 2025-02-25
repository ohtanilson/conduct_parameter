using Plots
using LinearAlgebra  # 線形代数の計算用
using Statistics     # 統計計算用

# Define the functions
P(Q, X) = X * exp(-Q)
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
plot_counter_example = plot(Q_range, Q -> P(Q, X1), label="P(Q_1,X_1)", 
     xlabel="Q", ylabel="Value", 
     title="P(Q,X) = X exp(-Q)")
plot!(Q_range, Q -> P(Q, X2), label="P(Q_2,X_2)")
hline!([P1], linestyle=:dash, color=:black, label="")
hline!([0], linestyle=:dash, color=:gray, label="", alpha=0.5)


# Plot ∂P/∂Q × Q for different X values
plot!(Q_range, Q -> ∂P∂Q(Q, X1)+ P(Q, X1), label="P(Q_1,X_1) + ∂P/∂Q×Q(Q,X_1)", linestyle=:dash)
plot!(Q_range, Q -> ∂P∂Q(Q, X2)+ P(Q, X2), label="P(Q_2,X_2) + ∂P/∂Q×Q(Q,X_2)", linestyle=:dash)

# Add points of interest
scatter!([Q1, Q2], [P1, P2], label="", color=:black)




# Save the plot
savefig(plot_counter_example, "../conduct_parameter/figuretable/plot_counter_example.png")

