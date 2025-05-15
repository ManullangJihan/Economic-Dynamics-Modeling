# Samuelson's Stability Diagram and Model Simulation

# Activate the Julia environment (adjust the path as needed)
include("../env/activate_env.jl")

using LinearAlgebra
using DifferentialEquations
using Distributions
using Plots
using Plots: plot, plot!

# Stability Diagram Section

# Define the range of k values (accelerator values) for plotting
k_end = 5
k1 = LinRange(0, k_end, 100)
k2 = LinRange(1, k_end, 100)
k3 = [0, k_end]

# Calculate the first stability boundary (b1) using the equation b1 = 4k / (1 + k)^2
b1 = 4 .* k1 ./ (1 .+ k1).^2

# Calculate the second stability boundary (b2) using the equation b2 = 1 / k
b2 = 1 ./ k2

# Define the third stability boundary (b3) as a constant b = 1
b3 = [1, 1]

# Plot the stability boundaries
plot(k1, b1, color=:black, linewidth=2, xlabel="k", ylabel="b",
xlim=(0, 5), ylim=(0, 1.2), label="k1")
plot!(k2, b2, color=:red, linewidth=2, label="k2")
plot!(k3, b3, color=:black, linewidth=2, linestyle=:dash, label="k3")

# Simulation of Samuelson's Model

# Set model parameters
b = 0.9  # Marginal propensity to consume
k = 0.2  # Accelerator (sensitivity of investment to output change)
G = 10   # Autonomous expenditure (constant external demand)

# Calculate the equilibrium value (fixed point)
Yfp = G / (1 - b)

# Initialize initial conditions with slight deviations from equilibrium
Y0 = Yfp - 0.1
Y1 = Yfp - 0.2

# Set the number of simulation steps
n = 8

# Initialize the output array (Y) with zeros
Y = zeros(n + 1)
Y[1] = Y0
Y[2] = Y1

# Simulate the difference equation of Samuelson's model
for t = 3:n+1
    Y[t] = b * (1 + k) * Y[t - 1] - b * k * Y[t - 2] + G
end

# Visualize the simulation result
plot(0:n, Y, color=:black, markersize=4, xlabel="t", ylabel="Yₜ", legend=false)
plot!([0, n], [Yfp, Yfp], linestyle=:dash, color=:black)
savefig("figures/Samuelson_model.png")