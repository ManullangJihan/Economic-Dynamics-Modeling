# Dornbusch Overshooting Model

# Activate environment and install dependencies defined in Project.toml
include("../env/activate_env.jl")

using LinearAlgebra               # For matrix operations
using DifferentialEquations       # To solve ODE systems numerically
using Plots                      # For plotting results

# Define the equilibrium exchange rate (constant)
const e_eq = 2.0

# Define the equilibrium price level (constant)
const p_eq = 2.0

# Define the system of differential equations describing dynamics around equilibrium
function eqn(du, u, p, t)
    # Unpack current state vector u into exchange rate e and price level p
    e, p = u

    # Define model parameters controlling system dynamics
    α = 1.0
    β = 1.0
    σ = 1.0
    θ = 0.5

    # Define the Jacobian matrix J based on parameters
    # This matrix defines how deviations from equilibrium evolve over time
    J = [0 1/α;                 # First row: change in e depends on p deviation
         θ*β -θ*(β + σ/α)]      # Second row: change in p depends on e and p deviations

    # Compute the time derivative vector dedp = J * (current deviation from equilibrium)
    dedp = J * [e - e_eq, p - p_eq]

    # Assign derivatives of e and p into du vector (output)
    du[1] = dedp[1]             # Rate of change of exchange rate e
    du[2] = dedp[2]             # Rate of change of price level p
end

# Set model parameters explicitly (matching those inside eqn function)
α = 1.0
β = 1.0
σ = 1.0
θ = 0.5

# Bundle parameters into a tuple (not used directly in eqn here, but good practice)
p = (α, β, σ, θ)

# Define simulation start and end times
t_start = 0.0
t_end = 5.0

# Create time span tuple for ODE solver
tspan = (t_start, t_end)

# Set initial conditions away from equilibrium to see dynamic response
e₀ = 3.0       # Initial exchange rate
p₀ = 1.0       # Initial price level
u₀ = [e₀, p₀]  # Initial state vector

# Compute Jacobian matrix explicitly (for reference or analysis)
J = [0 1/α; 
     θ*β -θ*(β + σ/α)
    ]

# Calculate initial derivative at starting point (not required for simulation, just example)
J * [e₀ - e_eq, p₀ - p_eq]

# Set solver tolerances to control numerical accuracy
reltol = 1e-6   # Relative tolerance
abstol = 1e-6   # Absolute tolerance

# Define the ODE problem with the system function, initial condition, and time span
problem = ODEProblem(eqn, u₀, tspan)

# Solve the ODE numerically using Tsit5 method with specified tolerances and save output at intervals of 0.1 time units
sol = solve(problem, Tsit5(), saveat=0.1, reltol=reltol, abstol=abstol)

# Extract time points from solution object
t = sol.t

# Convert solution object to array for easier plotting (rows: variables, columns: time steps)
sol = Array(sol)

# Plot exchange rate e(t) over time with thick black line
plot(t, sol[1, :], color=:black, linewidth=3.0, label="E",
     xlabel="Time")

# Add a red scatter point at initial exchange rate (time zero)
scatter!([0], [e₀], markercolor=:red, markersize=3, label=false)

# Plot price level p(t) over time with thick black line
plot(t, sol[2, :], color=:black, linewidth=3.0, label="P",
     xlabel="Time")

# Add a red scatter point at initial price level (time zero)
scatter!([0], [p₀], markercolor=:red, markersize=3, label=false)
savefig("figures/Dornbusch_model.png")