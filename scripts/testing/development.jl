# playing around with the odes 

using DifferentialEquations, ModelingToolkit, OrdinaryDiffEq, Plots, Latexify, ParameterizedFunctions
using ModelingToolkit: t_nounits as t, D_nounits as D

function lotka_volterra!(du, u, p, t)
    N₁, N₂ = u
    r₁, r₂, α₁₂, α₂₁, K₁, K₂ = p
    du[1] = dN₁ = r₁ * N₁ * (1 - ((N₁ + α₁₂ * N₂) / K₁))
    du[2] = dN₂ = r₂ * N₂ * (1 - ((N₂ + α₂₁ * N₁) / K₂))
end

function lotka_volterra_prop!(du, u, p, t)
    N₁, N₂ = u
    r₁, r₁ᵩ, r₂, α₁₂, α₂₁, K₁, K₂, ϕ = p
    dN₁ = ((ϕ * r₁ᵩ * N₁) * (1 - ϕ) * r₁ * N₁) * (1 - ((N₁ + α₁₂ * N₂) / K₁))
    dN₂ = r₂ * N₂ * (1 - ((N₂ + α₂₁ * N₁) / K₂))
end
# define also with parameterized functions
lv = @ode_def basicLV begin
    dN₁ = r₁ * N₁ * (1 - ((N₁ + α₁₂ * N₂) / K₁))
    dN₂ = r₂ * N₂ * (1 - ((N₂ + α₂₁ * N₁) / K₂))
end r₁ r₂ α₁₂ α₂₁ K₁ K₂
latexify(lv)

lv_prop = @ode_def lvProp begin
    dN₁ = ((ϕ * r₁ᵩ * N₁) * ((1 - ϕ) * r₁ * N₁)) * (1 - ((N₁ + α₁₂ * N₂) / K₁))
    dN₂ = r₂ * N₂ * (1 - ((N₂ + α₂₁ * N₁) / K₂))
end r₁ r₁ᵩ r₂ α₁₂ α₂₁ K₁ K₂ ϕ
latexify(lv_prop)


using DifferentialEquations, Optimization, OptimizationPolyalgorithms, SciMLSensitivity
using ForwardDiff, Plots

# Define experimental data
t_data = 0:10
x_data = [1.000 2.773 6.773 0.971 1.886 6.101 1.398 1.335 4.353 3.247 1.034]
y_data = [1.000 0.259 2.015 1.908 0.323 0.629 3.458 0.508 0.314 4.547 0.906]
xy_data = vcat(x_data, y_data)

# Plot the provided data
scatter(t_data, xy_data', label=["x Data" "y Data"])

# Setup the ODE function
function lotka_volterra!(du, u, p, t)
    x, y = u
    α, β, δ, γ = p
    du[1] = dx = α * x - β * x * y
    du[2] = dy = -δ * y + γ * x * y
end

# Initial condition
u0 = [1.0, 1.0]

# Simulation interval
tspan = (0.0, 10.0)


# LV equation parameter. p = [α, β, δ, γ]
pguess = [1.0, 1.2, 2.5, 1.2]

# Set up the ODE problem with our guessed parameter values
prob = ODEProblem(lotka_volterra!, u0, tspan, pguess)

# Solve the ODE problem with our guessed parameter values
initial_sol = solve(prob, saveat=1)

# View the guessed model solution
plt = plot(initial_sol, label=["x Prediction" "y Prediction"])
scatter!(plt, t_data, xy_data', label=["x Data" "y Data"])































ModelingToolkit.Latexify(lotka_volterra!)

# compute the Jacobian matrix
prob = ODEProblem(lotka_volterra!, [100, 100], (0, 1e5), (1, 1, 0.1, 0.1, 0.2, 0.2))
# modeltoolkit transformation 
@mtkbuild sys = modelingtoolkitize(prob)

prob_jac = ODEProblem(sys, [], (0.0, 1e5), jac=true)
lvjac = eval(generate_jacobian(sys))
plot(sol)