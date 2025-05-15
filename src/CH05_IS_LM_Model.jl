# IS-LM Model — Mathematical Explanation with Julia Implementation

# Activate environment and install dependencies defined in Project.toml
include("../env/activate_env.jl")

using LinearAlgebra
using DifferentialEquations
using Plots
using Plots: plot, plot!

# ------------------------
# PARAMETERS (Constants)
# ------------------------
cl = 1       # coefficient related to consumption/investment adjustment speed (c1)
c2 = 1       # coefficient related to money market adjustment speed (c2)
s = 0.5      # sensitivity of investment to output (Y) (marginal propensity to save)
a2 = 1       # sensitivity of investment to interest rate (r)
kl = 1       # sensitivity of money demand to output (k1)
k2 = 0.1     # sensitivity of money demand to interest rate (k2)
Ibar = 0.55  # autonomous investment (I_0)
Mbar = 1     # autonomous money supply (M_0)

# -----------------------------------------
# DEFINE COEFFICIENT MATRIX A AND VECTOR b
# System of linear ODEs: d/dt [Y, r] = A * [Y, r] + b
# -----------------------------------------
A = [-cl*s -cl*a2;     # First row: Effects on dY/dt from Y and r
     c2*kl -c2*k2]    # Second row: Effects on dr/dt from Y and r

b = [cl*Ibar,         # Constant term in dY/dt equation (investment term)
     -c2*Mbar]        # Constant term in dr/dt equation (money supply term)

# -----------------------------------------
# SIMULATION SETUP
# -----------------------------------------
t_start = 0
t_end = 200
tspan = (t_start, t_end)    # Time span for simulation

# Initial values for output (Y) and interest rate (r)
Y0 = 1
r0 = 0.01

# -----------------------------------------
# DEFINE SYSTEM OF DIFFERENTIAL EQUATIONS
# u is the state vector [Y, r], p is parameters (unused), t is time
# dy/dt = A*u + b
# -----------------------------------------
eqn(u, p, t) = A * u + b

# Define ODE problem
prob = ODEProblem(eqn, [Y0, r0], tspan)

# Solve the ODE numerically
sol = solve(prob, saveat=0.1) # Save solution at intervals of 0.1 time units

# Extract solution arrays
t = sol.t         # Time points
Z = Array(sol)    # Solution matrix with rows: [Y; r]

Y = Z[1, :]       # Output over time
r = Z[2, :]       # Interest rate over time

# -----------------------------------------
# CALCULATE EQUILIBRIUM POINT
# Equilibrium: d/dt [Y, r] = 0 => A * [Y_eq; r_eq] + b = 0
# Solve for equilibrium state: [Y_eq; r_eq] = -A^{-1} * b
# -----------------------------------------
Yr_eq = A \ (-b)   # Backslash operator solves linear system

Y_eq = Yr_eq[1]    # Equilibrium output
r_eq = Yr_eq[2]    # Equilibrium interest rate

# -----------------------------------------
# STABILITY ANALYSIS
# Calculate eigenvalues and eigenvectors of A
# Eigenvalues determine if system tends to equilibrium (stable) or diverges (unstable)
# -----------------------------------------
eigenvalues, eigenvecs = eigen(A)
lambda_1 = eigenvalues[1]
lambda_2 = eigenvalues[2]

println("Eigenvalues of A: ", lambda_1, ", ", lambda_2)

# -----------------------------------------
# PLOTTING RESULTS
# -----------------------------------------

# Plot output Y over time
plot(t, Y, color=:black, linewidth=3,
     xlim=(0, 20), ylim=(0.98, 1.04),
     xlabel="time (t)", ylabel="Output Y(t)",
     title="Output Y over time")

# Plot interest rate r over time
plot(t, r, color=:black, linewidth=3,
     xlim=(0, 20), ylim=(0.0, 0.08),
     xlabel="time (t)", ylabel="Interest Rate r(t)",
     title="Interest Rate r over time")

# Plot phase diagram: r vs Y showing system trajectory
plot(Y, r, color=:black, linewidth=3,
     xlim=(0.98,1.04), ylim=(0,0.08),
     xlabel="Output Y", ylabel="Interest Rate r",
     title="Phase Diagram of IS-LM Dynamics", label=false)

# Mark equilibrium point with a red dot
scatter!([Y_eq], [r_eq], color=:red, label="Equilibrium")

# Mark initial condition with a blue dot
scatter!([Y0], [r0], color=:blue, label="Initial State")

# Calculate and plot LM curve (Liquidity preference - Money supply equilibrium)
# LM curve formula from money market equilibrium: r = (k1/k2)*Y - Mbar/k2
Y_LM = [1, 1.008]                         # Select Y range for LM curve
r_LM = @. kl/k2*Y_LM - Mbar/k2           # Corresponding r values for LM curve

plot!(Y_LM, r_LM, color=:black, linestyle=:dash, label="LM Curve")

# Calculate and plot IS curve (Investment-Saving equilibrium)
# IS curve formula from goods market equilibrium: r = (Ibar/a2) - (s/a2)*Y
Y_IS = [0.98, 1.04]                      # Select Y range for IS curve
r_IS = @. Ibar/a2 - s/a2*Y_IS            # Corresponding r values for IS curve

plot!(Y_IS, r_IS, color=:black, linestyle=:dot, label="IS Curve")
savefig("figures/IS_LM_Dyanmics.png")