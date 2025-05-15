# Activate the Julia environment (adjust the path as needed)
include("../env/activate_env.jl")

# Load required libraries
using LinearAlgebra
using DifferentialEquations
using Plots
using Plots: plot, plot!


# ========================================
# Model Parameters
# ========================================

# s: Savings rate (proportion of output saved and reinvested)
const s = 0.8

# α: Output elasticity of capital (how output responds to capital changes)
const α = 0.5

# n: Population growth rate
const n = 0.5

# a: Technological growth rate
const a = 0.3

# δ: Depreciation rate of capital
const δ = 0.4


# ========================================
# Defining the Solow-Swan Differential Equation
# ========================================

# Differential equation for capital per worker (r)
function solow_system!(du, u, p, t)
    # Extract current value of r (capital per worker)
    r = u[1]

    # Calculate the rate of change of r (dr/dt)
    du[1] = drdt = s * r^α - (n + a + δ) * r
end


# ========================================
# Initial Condition and Time Setup
# ========================================

# Initial value of capital per worker (r)
r0 = 0.9
u0 = [r0]  # u0 is a vector since ODE solvers expect vectors

# Time range for simulation
t_start = 0  # Simulation start time
t_end = 10   # Simulation end time
tspan = (t_start, t_end)  # Define the time span


# ========================================
# Solving the Differential Equation
# ========================================

# Define the problem for the ODE solver
prob = ODEProblem(solow_system!, u0, tspan)

# Solve the problem using Tsit5() method (efficient ODE solver)
sol = solve(prob, Tsit5(), saveat=0.1)

# Extract time points and solution values
t = sol.t  # Time points
arr_sol = Array(sol)  # Convert solution to an array


# ========================================
# Calculating the Steady-State Value
# ========================================

# The equilibrium value of r (steady-state capital per worker)
r_eq = (s / (n + a + δ))^(1 / (1 - α))


# ========================================
# Plotting the Results
# ========================================

# Extract capital values (R) from the solution array
R = arr_sol[1, :]

# Plot the capital per worker over time
plot(t, R, color=:black, linewidth=3.0, label="R", xlabel="Time")

# Plot the equilibrium line
plot!([t_start, t_end], r_eq .* [1, 1], color=:black, linewidth=1, label=false)
savefig("figures/solow_swan_growth_model.png")