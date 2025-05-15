# Activate environment and install dependencies defined in Project.toml
include("../env/activate_env.jl")

# Load required libraries
using LinearAlgebra
using DifferentialEquations
using Plots
using Plots: plot, plot!

"""
Define Demand and Supply models
Demand depends on the current price (pₜ)
Dₜ = d₀ + d₁pₜ
Supply depends on the previous price (pₜ₋₁)
Sₜ = s₀ + s₁pₜ₋₁
"""

# Initial values for demand and supply parameters
# d0 is the base demand (demand when price is zero)
d0 = 4
# d1 is the price sensitivity of demand (negative means demand decreases with price)
d1 = -1.7

# s0 is the base supply (supply when price is zero)
s0 = 0.5
# s1 is the price sensitivity of supply (positive means supply increases with price)
s1 = 1.1

# Initial price and number of iterations for the simulation
p0 = 1.0    # Initial Price
n = 10      # Number of simulation steps

# Initialize price array with zeros and set the first value to initial price
p = zeros(Float64, n+1)
p[1] = p0

# Calculate price over time based on supply-demand balance
for t = 1:n
    # Price calculation using the supply-demand model
    p[t+1] = (s1/d1)*p[t] + (s0 - d0) / d1
end

# Plot the price adjustment over time
plot(0:n, p, color=:black, linewidth=2, legend=false)
scatter!(0:n, p, color=:red)

# Change the supply sensitivity to see its effect
s1 = 1.9
p[1] = p0

# Re-calculate price with new supply sensitivity
for t = 1:n
    p[t+1] = (s1/d1)*p[t] + (s0 - d0) / d1
end


# Plot the adjusted scenario
plot(0:n, p, color=:black, linewidth=2, label="Price(t)", xlabel="t", ylabel="Price")
scatter!(0:n, p, color=:red, label=false)

# Plotting the demand and supply lines
plot([0, d0], [-d0/d1, 0], color=:black, label="Demand", xlabel="Quantity (q)", ylabel="Price (p)")
plot!([s0, d0], [0, (d0-s0)/s1], color=:black, label="Supply")

# Plotting the equilibrium point (intersection of demand and supply)
scatter!(
    [(d1*s0-d0*s1)/(d1-s1)], [(s0-d0)/(d1-s1)], color=:red,
    markersize=6, label="Equilibrium"
)

savefig("figures/demand_and_supply_equilibrium.png")