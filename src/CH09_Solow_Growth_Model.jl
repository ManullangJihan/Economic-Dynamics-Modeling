# Solow Growth Model

# Load environment and dependencies
include("../env/activate_env.jl") # Activate the Julia environment

# Import required libraries
using LinearAlgebra
using Plots
using Plots: plot, plot!

# Model Parameters
s = 0.8   # Savings rate (fraction of output saved and invested)
α = 0.5   # Capital share in the production function
n = 0.5   # Population growth rate (labor growth)
a = 0.3   # Technological growth rate
δ = 0.4   # Depreciation rate of capital

# Define the capital per effective worker (r) range
r1 = collect(LinRange(0,1,100)) # A range of r values from 0 to 1 for smooth plotting
r2 = [0, 1]                     # Simple 2-point range for straight line plotting

# Calculate the steady-state value of capital per effective worker
r_eq = (s / (n + a + δ))^(1 / (1 - α)) # Steady-state capital value calculation

# Plotting the Solow Model Graph
plot(r1, s .* r1.^α, color=:black, linewidth=3.0, label="Investment (y1)", xlabel="r") # Investment curve
plot!(r2, (n + a + δ) .* r2, color=:black, linewidth=3, label="Depreciation (y2)") # Depreciation line

# Dashed lines for visualizing the steady-state point
plot!(r_eq .* [1, 1], s * r_eq^(α) .* [0, 1], color=:black, linestyle=:dash, label=false) # Vertical line at r_eq
plot!(r_eq .* [0, 1], s * r_eq^(α) .* [1, 1], color=:black, linestyle=:dash, label=false) # Horizontal line at steady-state output
savefig("figures/solow_growth_model.png")