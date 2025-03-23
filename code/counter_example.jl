using Plots
using LinearAlgebra  # 線形代数の計算用
using Statistics     # 統計計算用

## Define the functions
θ = 0
θ_alternative = 1

P(Q, X) = X * exp(-Q) + 1
∂P∂Q(Q, X) = - X * exp(-Q)
MR(Q, X, θ) = P(Q, X) + θ * ∂P∂Q(Q, X) * Q


G1 = MR(Q1, X1, θ)
G2 = MR(Q2, X2, θ)

# Test points
Q1, X1 = 1.0, 2.0
Q2, X2 = 2.0, 2exp(1)  # ≈ 5.437

## Visualization

Q_range = 0:0.1:5
X_values = [2.0, 2exp(1)]

# Plot P(Q,X) for different X values
plot_counter_example = plot(label="MR(Q_1,X_1;0)", xlabel="Q", ylabel="Value", 
     title="P(Q,X) = X exp(-Q) + 1", color = :blue)
plot!(Q_range, Q -> MR(Q, X1, θ), label="MR(Q_1,X_1;0)", color = :red)
plot!(Q_range, Q -> MR(Q, X2, θ), label="MR(Q_2,X_2;0)", color = :blue)
hline!([0], linestyle=:dash, color=:gray, label="", alpha=0.5)
hline!([G1], linestyle=:dash, color=:black, label="")
# Plot ∂P/∂Q × Q for different X values
plot!(Q_range, Q -> MR(Q, X1, θ_alternative), label="MR(Q_1,X_1;1)", linestyle=:dash, color = :blue)
plot!(Q_range, Q -> MR(Q, X2, θ_alternative), label="MR(Q_2,X_2;1)", linestyle=:dash, color = :red)

# Add points of interest
scatter!([Q1, Q2], [G1, G2], label="", color=:black)

display(plot_counter_example)

# Save the plot
savefig(plot_counter_example, "../conduct_parameter/figuretable/plot_counter_example.pdf")






########################################################

## Specific function in Theorem 

########################################################

θ = 1/2
θ_alternative = 1/3

# Define the functions
P(Q, X) = Q^(-1/θ) * X 
∂P∂Q(Q, X) = (-1/θ) * Q^(-1/θ-1) * X 

MR(Q, X, θ) = P(Q, X) + θ * ∂P∂Q(Q, X) * Q


Q_range = 0.5:0.01:1
X_values = [0.1, 0.11]

# Plot P(Q,X) for different X values
plot_counter_example = plot(xlabel="Q", ylabel="Value", title="P(Q,X) = Q^(-1/θ) * X", color = :blue)
plot!(Q_range, Q -> MR(Q, X1, θ), label="MR(Q_1,X_1;1/2)", color = :blue)
plot!(Q_range, Q -> MR(Q, X2, θ), label="MR(Q_2,X_2;1/2)", color = :red)
hline!([0], linestyle=:dash, color=:gray, label="", alpha=0.5)

# Plot ∂P/∂Q × Q for different X values
plot!(Q_range, Q -> MR(Q, X1, θ_alternative), label="MR(Q_1,X_1;1/3)", linestyle=:dash, color = :blue)
plot!(Q_range, Q -> MR(Q, X2, θ_alternative), label="MR(Q_2,X_2;1/3)", linestyle=:dash, color = :red)


display(plot_counter_example)

# Save the plot
savefig(plot_counter_example, "../conduct_parameter/figuretable/plot_counter_example_specific_function.pdf")



########################################################

## Specific function in Theorem 

########################################################

θ = 1/2
θ_alternative = 1/3

# Define the functions
P(Q, X) = X * exp(-Q)+1
∂P∂Q(Q, X) = - X * exp(-Q)  # This is ∂P/∂Q × Q

MR(Q, X, θ) = P(Q, X) + θ * ∂P∂Q(Q, X) * Q

Q1, X1 = 1.0, 2.0
Q2, X2 = 2.0, 2exp(1)  # ≈ 5.437

G1 = MR(Q1, X1, θ)
G2 = MR(Q2, X2, θ)

Q_range = 0:0.01:3
X_values = [0.1, 0.11]

# Plot P(Q,X) for different X values
plot_counter_example = plot(xlabel="Q", ylabel="Value", title="P(Q,X) = X exp(-Q) + 1", color = :blue)
plot!(Q_range, Q -> MR(Q, X1, θ), label="MR(Q_1,X_1;1/2)", color = :blue)
plot!(Q_range, Q -> MR(Q, X2, θ), label="MR(Q_2,X_2;1/2)", color = :red)
hline!([0], linestyle=:dash, color=:gray, label="", alpha=0.5)
hline!([G1], linestyle=:dash, color=:black, label="")
plot!(Q_range, Q -> MR(Q, X1, θ_alternative), label="MR(Q_1,X_1;8/9)", linestyle=:dash, color = :blue)
plot!(Q_range, Q -> MR(Q, X2, θ_alternative), label="MR(Q_2,X_2;8/9)", linestyle=:dash, color = :red)


display(plot_counter_example)



########################################################

## Specific function in Theorem 

########################################################



θ = 1/2
θ_alternative = 1/3

# Define the functions
P(Q, X) = Q^(-1/4) * X 
∂P∂Q(Q, X) = (-1/4) * Q^(-5/4) * X 

MR(Q, X, θ) = P(Q, X) + θ * ∂P∂Q(Q, X) * Q


Q_range = 0.1:0.01:1
X_values = [0.1, 0.1]

# Plot P(Q,X) for different X values
plot_counter_example = plot(xlabel="Q", ylabel="Value", title="P(Q,X) = Q^(-1/θ) * X", color = :blue)
plot!(Q_range, Q -> MR(Q, X1, θ), label="MR(Q_1,X_1;1/2)", color = :blue)
plot!(Q_range, Q -> MR(Q, X2, θ), label="MR(Q_2,X_2;1/2)", color = :red)
hline!([0], linestyle=:dash, color=:gray, label="", alpha=0.5)

# Plot ∂P/∂Q × Q for different X values
plot!(Q_range, Q -> MR(Q, X1, θ_alternative), label="MR(Q_1,X_1;1/3)", linestyle=:dash, color = :blue)
plot!(Q_range, Q -> MR(Q, X2, θ_alternative), label="MR(Q_2,X_2;1/3)", linestyle=:dash, color = :red)


display(plot_counter_example)




########################################################

## The linear demand model 

########################################################


θ = 0
θ_alt = 1

α_0, α_1, α_2 = 5, 2, 1
β_0_alt, β_1_alt, β_2_alt = 1, 1,1
β_0, β_1, β_2 = β_0_alt, α_1 + β_1_alt, β_2_alt

α = [α_0, α_1, α_2]
β = [β_0, β_1, β_2]
β_alt = [β_0_alt, β_1_alt, β_2_alt]

P(Q, X_d) = α_0 - α_1 * Q + α_2 * X_d
∂P∂Q(Q, X_d) = -α_1

MC(Q,X_s) = β_0 + β_1 * Q + β_2 * X_s
MC_alt(Q,X_s) = β_0_alt + β_1_alt * Q + β_2_alt * X_s

Q_equ(X_d, X_s, θ, α, β) = (α[1] - β[1] + α[3]X_d - β[3]X_s) / ((θ +1)*α[2] + β[2])

MR(Q, X, θ) = P(Q, X) + θ * ∂P∂Q(Q, X) * Q

P_supply(Q, X_s, X_d, θ) = MC(Q, X_s) - θ * ∂P∂Q(Q, X_d) * Q

F(Q, X_d, X_s) = MR(Q, X_d, θ) - MC(Q, X_s)
F_alt(Q, X_d, X_s) = MR(Q, X_d, θ_alt) - MC_alt(Q, X_s)



# MR(Q1,X1_d,0) = MR(Q2,X2_d,0)となるようなX1_d, X2_d, X1_s, X2_sを求める
X1_d, X1_s = 1, 1
X2_d, X2_s = 2, (-(-(α[2] * α[3])/(α[2] + β[2])  + α[3]) +(α[2] * β[3])/(α[2] + β[2])) * (α[2] + β[2])/(α[2] * β[3])


Q1 = Q_equ(X1_d, X1_s, θ, α, β)
Q2 = Q_equ(X2_d, X2_s, θ, α, β)

P1 = P(Q1, X1_d)
P2 = P(Q2, X2_d)


G1 = MR(Q1, X1_d, θ)
G2 = MR(Q2, X2_d, θ)

G1_alt = MR(Q1, X1_d, θ_alt)
G2_alt = MR(Q2, X2_d, θ_alt)

P_supply1_true = P_supply(Q1, X1_s, X1_d, θ)
P_supply2_true = P_supply(Q2, X2_s, X2_d, θ)

P_supply1_alt = P_supply(Q1, X1_s, X1_d, θ_alt)
P_supply2_alt = P_supply(Q2, X2_s, X2_d, θ_alt)

F1 = F(Q1, X1_d, X1_s)
F2 = F(Q2, X2_d, X2_s)

F1_alt = F_alt(Q1, X1_d, X1_s)
F2_alt = F_alt(Q2, X2_d, X2_s)

Q_range = 0:0.01:1.5
X_values = [0.1, 0.1]


plot_market_equilibrium = plot(xlabel="Q", ylabel="Value", title="Market Equilibrium", legend=:outertopright)
plot!(Q_range, Q -> P(Q, X1_d), label="P(Q1,X1_d)", color = :red)
plot!(Q_range, Q -> MC(Q, X1_s), label="MC(Q,X_s)", color = :blue)
# for θ_alt
plot!(Q_range, Q -> MR(Q, X1_d, θ_alt), label="MR(Q1,X1_d;1)", color = :red)
plot!(Q_range, Q -> MC_alt(Q, X1_s), label="MC_alt(Q,X_s)", color = :blue)
scatter!([Q1], [P1], label="", color=:black)
vline!([Q1], linestyle=:dash, color=:gray, label="")

plot!(Q_range, Q -> P(Q, X2_d), label="P(Q2,X2_d)", color = :red, linestyle=:dash)
plot!(Q_range, Q -> MC(Q, X2_s), label="MC(Q,X_s)", color = :blue, linestyle=:dash)
# for θ_alt
plot!(Q_range, Q -> MR(Q, X2_d, θ_alt), label="MR(Q2,X2_d;1)", color = :red, linestyle=:dash)
plot!(Q_range, Q -> MC_alt(Q, X2_s), label="MC_alt(Q,X_s)", color = :blue, linestyle=:dash)
scatter!([Q2], [P2], label="", color=:black)
vline!([Q2], linestyle=:dash, color=:gray, label="")
hline!([G1], linestyle=:dash, color=:gray, label="")
scatter!([Q1, Q2], [G1_alt, G2_alt], label="", color=:black, 
        markerstrokecolor=:black,    # 枠の色
        markercolor=:white)          # 内部を白に
hline!([G1_alt, G2_alt], linestyle=:dash, color=:gray, label="")





plot!(Q_range, Q -> F(Q, X1_d, X1_s), label="F(Q,X_d,X_s)", color = :red)
plot!(Q_range, Q -> F_alt(Q, X1_d, X1_s), label="F_alt(Q,X_d,X_s)", color = :blue, linestyle=:dash)
scatter!([Q1], [P1], label="", color=:black)
vline!([Q1], linestyle=:dash, color=:black, label="")

plot!(Q_range, Q -> F(Q, X2_d, X2_s), label="F(Q,X_d,X_s)", color = :red)
plot!(Q_range, Q -> F_alt(Q, X2_d, X2_s), label="F_alt(Q,X_d,X_s)", color = :blue, linestyle=:dash)
scatter!([Q2], [P2], label="", color=:black)
vline!([Q2], linestyle=:dash, color=:black, label="")
hline!([0], linestyle=:dash, color=:gray, label="", alpha=0.5)


plot_counter_example = plot(xlabel="Q", ylabel="Value", title="P(Q,X) = Q^(-1/θ) * X")
plot!(Q_range, Q -> P(Q, X1_d), label="P(Q,X_d;1/2)", color = :blue)
plot!(Q_range, Q -> MC(Q, X1_s), label="MC(Q,X_s;1/2)", color = :blue)
plot!(Q_range, Q -> MR(Q, X1_d, θ), label="MR(Q,X_d;1/2)", color = :red)




