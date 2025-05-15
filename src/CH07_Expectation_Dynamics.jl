# Activate environment and install dependencies defined in Project.toml
include("../env/activate_env.jl")

# Import necessary packages
using LinearAlgebra          
using Distributions
using Plots     
using Plots: plot, plot!

# -------------------------
# 1. Naive Expectation simulation
# -------------------------

α = 0.5   # Parameter alpha (affects dynamics)
γ = 0.5   # Parameter gamma (affects dynamics)
m = 10    # Constant term (target/fixed point)
p0 = 11   # Initial value of the sequence p at time 0
n = 10    # Number of time steps to simulate

p = zeros(n+1)   # Preallocate array to store p values from t=0 to t=n
p[1] = p0        # Set initial condition p(0) = p0

# Iterate over time steps to compute p(t) for t=1,...,n
for k = 1:n
    # Update rule derived from difference equation:
    # p(k+1) = (-α*γ)/(1 - α*γ) * p(k) + m/(1 - α*γ)
    p[k+1] = -α*γ / (1 - α*γ) * p[k] + m / (1 - α*γ)
end

# Plot the sequence p over time
plot(0:n, p, color=:black, linewidth=3.0, legend=false,
     xlabel="t", ylabel="p(t)")
scatter!(0:n, p, color=:green, markersize=5.0)

# -------------------------
# 2. Perfect Foresight simulation
# -------------------------

α = 4.0    # New alpha parameter (larger, changes dynamics)
m = 10     # Constant term remains the same
p0 = 11    # Initial value of p at time 0
n = 10     # Number of time steps

p = zeros(n+1)  # Preallocate p array
p[1] = p0       # Initial condition p(0) = p0

# Iterate to compute p(t) using perfect foresight update rule:
# p(k+1) = (1 + α)/α * p(k) - m/α
for k = 1:n
    p[k+1] = (1 + α) / α * p[k] - m / α
end

# Plot the sequence p(t) over time
plot(0:n, p, color=:black, linewidth=3.0, legend=false,
     xlabel="t", ylabel="p(t)")
scatter!(0:n, p, color=:green, markersize=5.0)

# -------------------------
# 3. Rational Expectation Simulation with noise
# -------------------------

α = 4.0    # Alpha parameter
γ = 0.5    # Gamma parameter (not used explicitly here)
m = 10     # Constant term (target/fixed point)
p0 = 10    # Initial condition for p(0)
n = 10     # Number of time steps

p = zeros(n+1)  # Preallocate p array
p[1] = p0       # Set initial value p(0) = p0

d = Normal(1, 1)  # Define Normal distribution with mean=1, std=1 for noise e(k)

# For each time step, sample noise and compute p(k+1)
for k = 1:n
    e = rand(d)    # Draw random noise sample e_k ~ N(1,1)
    # Update rule with noise:
    # p(k+1) = (1/(1+α)) * e_k + m
    p[k+1] = 1/(1 + α) * e + m
end

pfp = copy(m)  # Fixed point/reference level for plotting

# Plot noisy sequence p(t) with fixed point dashed line
plot(0:n, p, color=:black, linewidth=3.0, legend=false,
     xlabel="t", ylabel="p(t)")
scatter!(0:n, p, color=:green, markersize=5.0)
plot!([0, n], [pfp, pfp], linestyle=:dash, linewidth=2.0)
savefig("figures/rational_expectation_with_noise.png")