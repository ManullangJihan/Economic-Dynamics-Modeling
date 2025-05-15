# Activate environment and install dependencies defined in Project.toml
include("../env/activate_env.jl")

using LinearAlgebra
using DifferentialEquations
using Plots
using Plots: plot, plot!

# Parameters
b = 1                   # primary deficit to GDP ratio
r = 2                   # real interest rate
g = 3                   # growth rate of real GDP
p = [b, r, g]           

# Time span and initial condition
t_start = 0.0
t_end = 10.0
d0 = -1.0                # initial debt-to-GDP ratio
tspan = (t_start, t_end)

# ODE function: d'(t) = b + (r - g) * d(t)
function debt_eq!(du, u, p, t)
    b, r, g = p
    d = u[1]
    du[1] = b + (r - g) * d
end

# Define and solve ODE problem
prob = ODEProblem(debt_eq!, [d0], tspan, p)
sol = solve(prob, saveat=0.1, reltol=1e-6, abstol=1e-6)

# Extract solution and time points
d = Array(sol)[:]
t = sol.t

# Calculate equilibrium debt ratio d_eq = b / (g - r)
d_eq = b / (g - r)

# Plot results
plot(t, d, color=:black, linewidth=2, xlabel="Time (t)", ylabel="Debt-to-GDP Ratio d(t)", legend=false)
plot!([t_start, t_end], [d_eq, d_eq], linestyle=:dash, color=:red)

# Display equilibrium value on the plot title
title!("Debt Ratio over Time (Equilibrium d_eq = $(round(d_eq, digits=3)))")
savefig("figures/debt_to_gdp_ratio.png")