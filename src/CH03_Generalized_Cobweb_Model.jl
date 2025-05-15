# Generalized Cobweb Model
# Demand and supply with normal price expectation rule

# Activate environment and install dependencies defined in Project.toml
include("../env/activate_env.jl")

# Importing necessary libraries
using LinearAlgebra
using DifferentialEquations
using Plots
using Plots: plot, plot!

# Defining demand and supply parameters
# Demand coefficients
d0 = 4    # Intercept of demand
d1 = -1.7 # Slope of demand (negative for downward sloping)

# Supply coefficients
s0 = 0.5  # Intercept of supply
s1 = 1.9  # Slope of supply (positive for upward sloping)

# Adjustment factor (0 to 1) for price expectations
c = 0.8

# Initial price
p0 = 1

# Calculating equilibrium price (where demand equals supply)
pN = (s0 - d0) / (d1 - s1)

# Setting the number of iterations (time steps)
n = 10

# Initializing price array for storing simulated prices
p = zeros(Float64, n + 1)
p[1] = p0  # Initial price

# Simple Cobweb Model Simulation (without price expectations)
for k in 1:n
    # Calculating next price using the cobweb formula
    p[k+1] = (s1 / d1) * (1 - c) * p[k] + (s1 * c * pN + s0 - d0) / d1
end

# Plotting the simple cobweb model result
plot(0:n, p, color=:black, linewidth=3, legend=false)
scatter!(0:n, p, color=:grey, xlabel="t", ylabel="p(t)")


# Cobweb Model with Price Expectation (Adaptive)
# Re-initializing price arrays
p = zeros(n + 1)
pk = zeros(n + 1)
pe = zeros(n + 1)

# Initial price setup
p[1] = p0
pk[1] = p0

# Adaptive Expectations Simulation
for k in 1:n
    # Calculating expected price based on previous price and equilibrium
    pe[k+1] = p[k] + c * (pN - p[k])
    
    # Calculating next price using expected price
    p[k+1] = (s1 / d1) * pe[k+1] + (s0 - d0) / d1
end

# Calculating Demand and Supply for visualization
D = d0 .+ d1 .* p      # Demand values based on prices
S = s0 .+ s1 .* pe     # Supply values based on expected prices

# Setting up cobweb plot
p_obj = plot(
    [0, d0], [-d0/d1, 0], color=:black, legend=false,
    xlim=(2.28, 2.38),
    ylim=(0.94, 1.01),
    xlabel="Quantity (q)", ylabel="Price (p)"
)

# Plotting demand and supply curves
plot!(p_obj, [s0, d0], [0, (d0 - s0 / s1)], color=:black)

# Drawing the cobweb diagram (connecting points)
for k in 2:n
    # Connecting demand and supply points
    plot!(p_obj, [D[k-1], S[k]], [p[k-1], pe[k]], color=:black)
    scatter!(p_obj, [D[k-1], S[k]], [p[k-1], pe[k]], color=:black)

    # Alternating between dashed lines for the cobweb
    if k % 2 == 0 && k > 1
        plot!(p_obj, [S[k], D[k]], [pe[k], p[k]], color=:black, linestyle=:dash)
    elseif k > 1
        plot!(p_obj, [D[k], S[k]], [p[k], pe[k]], color=:black, linestyle=:dash)
    end
end

# Displaying the final cobweb diagram
display(p_obj)
savefig("figures/cobweb_diagram.png")