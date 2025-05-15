# Activate environment and install dependencies defined in Project.toml
include("../env/activate_env.jl")

# Importing required libraries
using LinearAlgebra
using DifferentialEquations
using Plots
using Plots: plot, plot!

# Defining model parameters
# 'a' is the autonomous consumption (fixed amount of consumption when income is zero)
a = 10

# 'b' is the marginal propensity to consume (MPC), which is the fraction of income spent on consumption
b = 0.8

# 'Ibar' is the autonomous investment (fixed amount of investment)
Ibar = 20

# 'Gbar' is the autonomous government spending (fixed amount)
Gbar = 20

# Setting the initial condition
# 'Y0' is the initial national income value
Y0 = 300  # Initial value (this line is redundant)
Y0 = 10   # Correct initial value

# Calculating the fixed point (equilibrium income)
# This is derived from the formula Y_fp = (a + Ibar + Gbar) / (1 - b)
Yfp = (a + Ibar + Gbar) / (1 - b)

# Setting simulation parameters
n = 30  # Number of iterations for the simulation

# Initializing a vector to store income values over time
Y = zeros(n+1)  # Vector of zeros, size n+1 to include initial value
Y[1] = Y0       # Setting the initial value

# Simulation loop (iteratively calculating national income)
for m = 1:n
    # Calculating the next value using the Keynesian model equation
    Y[m+1] = b * Y[m] + a + Ibar + Gbar
end

# Plotting the results
plot(0:n, Y, color=:black, xlabel="t", ylabel="Y(t)", legend=false)

# Adding scatter points for better visualization
scatter!(0:n, Y, color=:black)

# Plotting the equilibrium (fixed point) line
plot!([0,n], [Yfp, Yfp], color=:black, linestyle=:dash)
savefig("figures/keynesian_multiplier.png")

# The commented section below is an analytical solution using a geometric series approach
# Y = zeros(n+1)
# Y[1] = Y0
# for m = 1:n
#     Y[m+1] = (Y0 - Yfp) * b^m + Yfp
# end