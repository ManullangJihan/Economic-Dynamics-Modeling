# Endogenous Growth Model

# Load environment and dependencies
include("../env/activate_env.jl") # Activate the Julia environment

using LinearAlgebra # For matrix and vector calculations
using DifferentialEquations # For solving differential equations
using Distributions # For probability distributions (not used here but loaded)
using Plots # For visualization
using Plots: plot, plot! # Specific functions from Plots

# Parameters (Economic Model Parameters)
const α = 0.5 # Capital elasticity (0 < α < 1)
const β = 0.2 # Influence of capital on knowledge growth
const γ = 2.0 # Base growth rate of knowledge
const θ = 0.5 # Self-reinforcement of knowledge
const n = 1.0 # Constant external growth factor
p = (α, β, γ, θ, n) # Collect parameters into a tuple (optional)

# Equilibrium Points Calculation
# These are the theoretical solutions where the system can stabilize
eq1 = [0, 0] # Trivial equilibrium
eq2 = [n, 0] # Only capital grows
eq3 = [0, γ*n/(1-θ)] # Only knowledge grows
# Both capital and knowledge grow at a balanced rate
eq4 = [n*(β+γ)/(1-(β+θ))+n, n*(β+γ)/(1-(β+θ))]

# Jacobian Matrix Calculation (Stability Analysis)
J4 = [
    (1-α)*(eq4[2]+2-2*eq4[1]) (1-α)*eq4[1];
    β*eq4[2] β*eq4[1]+γ*n+2*(θ-1)*eq4[2]
]

# Initial Conditions for the ODE System
g_K0 = 14.0 # Initial value for capital growth
g_A0 = 1.0  # Initial value for knowledge growth
y0 = [g_K0, g_A0] # Combined as a vector for the solver

# Function that Defines the ODE System
function endogenous_system(du, u, p, t)
    g_K, g_A = u # Unpack the current values
    # Differential equations for capital and knowledge
    du[1] = (1-α) * (g_A + n - g_K) * g_K
    du[2] = (β * g_K + γ * n + (θ - 1) * g_A) * g_A
end

# Simulation Time Setup
t_start = 0.0 # Start time
t_end = 10.0 # End time
tspan = (t_start, t_end) # Time span for the simulation

# Solver Setup
reltol = 1e-6 # Relative tolerance for high accuracy
abstol = 1e-6 # Absolute tolerance for high accuracy

# Setting Up the ODE Problem
problem = ODEProblem(endogenous_system, y0, tspan)

# Solving the ODE using Tsitouras 5/4 Runge-Kutta method
sol = solve(problem, Tsit5(), saveat=0.1, reltol=reltol, abstol=abstol)

# Extracting Solutions
t = sol.t # Time points
g_K = sol[1, :] # Capital growth over time
g_A = sol[2, :] # Knowledge growth over time

# Eigenvalues and Eigenvectors for Stability Analysis
eigenvals, eigenvecs = eigen(J4)
mu4_1 = eigenvals[1] # First eigenvalue
mu4_2 = eigenvals[2] # Second eigenvalue

# Visualization (Figure 1)
p1 = plot(t, g_K, color=:black, linewidth=2, label="g_K(t)")
plot!(p1, t, g_A, color=:black, linewidth=2, label="g_A(t)")

# Phase Plot (Figure 2)
p2 = plot(g_K, g_A, color=:black, linewidth=2, label=false, xlabel="g_K(t)", ylabel="g_A(t)")
scatter!(p2, [eq1[1]], [eq1[2]], markercolor=:black, label="Equilibrium 1")
scatter!(p2, [eq2[1]], [eq2[2]], markercolor=:red, label="Equilibrium 2")
scatter!(p2, [eq3[1]], [eq3[2]], markercolor=:green, label="Equilibrium 3")
scatter!(p2, [eq4[1]], [eq4[2]], markercolor=:yellow, label="Equilibrium 4")
savefig("figures/Endogenous_Growth_Model.png")