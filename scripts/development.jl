# playing around with the odes 

using DifferentialEquations, ModelingToolkit, OrdinaryDiffEq, Plots
using ModelingToolkit: t_nounits as t, D_nounits as D

function lotka_volterra!(du, u, p, t)
    N₁, N₂ = u
    r₁, r₂, α₁₁, α₁₂, α₂₁, α₂₂ = p
    du[1] = dN₁ = N₁ * (r₁ - (α₁₁ * N₁) - (α₁₂ * N₂))
    du[2] = dN₂ = N₂ * (r₂ - (α₂₁ * N₁) - (α₂₂ * N₂))
end

# compute the Jacobian matrix
prob = ODEProblem(lotka_volterra!, [100, 100], (0, 1e5), (1, 1, 0.1, 0.1, 0.2, 0.2))
# modeltoolkit transformation 
@mtkbuild sys = modelingtoolkitize(prob)

prob_jac = ODEProblem(sys, [], (0.0, 1e5), jac=true)
lvjac = eval(generate_jacobian(sys))
plot(sol)