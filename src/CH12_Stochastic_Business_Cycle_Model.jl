# Stochastic Business Cycle Model

# Activate the Julia environment (adjust the path as needed)
include("../env/activate_env.jl")

using LinearAlgebra
using SparseArrays
using Plots

# Model Parameters (single set)
α = 0.6        # Capital elasticity of output
β = 0.975      # Discount factor
δ = 0.025      # Depreciation rate
σ = 0.5        # Risk aversion / curvature parameter
A = 1          # Technology parameter

# Initial capital stock
a₀ = 10

# Time horizon
T = 100
n = 0:1:T+1
N = length(n)

# Initial guess for capital path (linear decreasing)
k₀ = @. a₀ * (1 - n/(T+1))

# Function to compute residual vector G(k)
function gfun(k, a0, A, alpha, beta, delta, sigma)
    # Slice k into overlapping windows
    k1_temp = k[1:end-2]
    k2_temp = k[2:end-1]
    k3_temp = k[3:end]

    # Residuals for Euler equation system
    R1 = @. A * k1_temp^alpha + (1 - delta) * k1_temp - k2_temp
    R2 = @. A * k2_temp^alpha + (1 - delta) * k2_temp - k3_temp

    # Marginal product derivative term
    Q = @. 1 + alpha * A * k2_temp^(alpha - 1) - delta

    # Boundary conditions
    FT0 = k[1] - a0    # initial capital constraint
    FT1 = k[end]       # terminal condition (can be zero or free)

    # Euler residual vector
    F = @. R1^(-sigma) - beta * Q * R2^(-sigma)

    # Compose full residual vector: first BC, Euler eqs, last BC
    G = vcat(FT0, F, FT1)
    return G
end

# Function to compute Jacobian matrix dG/dk
function dGfun(k, A, alpha, beta, delta, sigma)
    k1_temp = k[1:end-2]
    k2_temp = k[2:end-1]
    k3_temp = k[3:end]

    R1 = @. A * k1_temp^alpha + (1 - delta) * k1_temp - k2_temp
    R2 = @. A * k2_temp^alpha + (1 - delta) * k2_temp - k3_temp

    dR1_dt = @. alpha * A * k1_temp^(alpha - 1) + 1 - delta
    dR1_dt1 = -1
    dR2_dt1 = @. alpha * A * k2_temp^(alpha - 1) + 1 - delta
    dR2_dt2 = -1

    Q = @. 1 + alpha * A * k2_temp^(alpha - 1) - delta
    dQ_t1 = @. alpha * (alpha - 1) * A * k2_temp^(alpha - 2)

    dG_t = @. -sigma * R1^(-sigma - 1) * dR1_dt
    dG_t1 = @. -sigma * R1^(-sigma - 1) * dR1_dt1 - beta * dQ_t1 * R2^(-sigma) - beta * Q * (-sigma) * R2^(-sigma - 1) * dR2_dt1
    dG_t2 = @. -beta * Q * (-sigma) * R2^(-sigma - 1) * dR2_dt2

    N = length(k)

    # Initialize Jacobian as sparse matrix N x N
    dGmat = spzeros(N, N)

    # Boundary conditions derivatives (first and last rows)
    dGmat[1, 1] = 1
    dGmat[end, end] = 1

    # Fill in the middle rows corresponding to Euler equations
    for i in 1:(N-2)
        dGmat[i+1, i] = dG_t[i]      # diagonal below main (offset -1)
        dGmat[i+1, i+1] = dG_t1[i]  # main diagonal
        if i < N-2
            dGmat[i+1, i+2] = dG_t2[i]  # diagonal above main (offset +1)
        end
    end

    return dGmat
end

# Newton-Raphson iteration to solve for k
Kold = k₀
Knew = copy(k₀) 
err = 1.0

while err > 1e-4
    G = gfun(Kold, a₀, A, α, β, δ, σ)
    dG = dGfun(Kold, A, α, β, δ, σ)
    # Newton step: Knew = Kold - Jacobian \ Residuals
    Knew = Kold - dG \ G
    err = norm(Knew - Kold)
    Kold = Knew
    println("Error: ", err)
end

# Plot the solution capital path
plot(n, Knew, label="Capital Path", xlabel="Time", ylabel="Capital k", linewidth=2, legend=false)
savefig("figures/stochastic_business_cycle_model.png")

