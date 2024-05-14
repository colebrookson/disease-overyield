# playing around with the odes 

using DifferentialEquations, ModelingToolkit, OrdinaryDiffEq, Plots, Latexify,
    ParameterizedFunctions
using ModelingToolkit: t_nounits as t, D_nounits as D

function lotka_volterra!(du, u, p, t)
    N₁, N₂ = u
    r₁, r₂, α₁₂, α₂₁ = p
    du[1] = dN₁ = r₁ * N₁ * (1 - ((N₁ + α₁₂ * N₂) / K₁))
    du[2] = dN₂ = r₂ * N₂ * (1 - ((N₂ + α₂₁ * N₁) / K₂))
end

# define also with parameterized functions



ModelingToolkit.Latexify(lotka_volterra!)

# compute the Jacobian matrix
prob = ODEProblem(lotka_volterra!, [100, 100], (0, 1e5), (1, 1, 0.1, 0.1, 0.2, 0.2))
# modeltoolkit transformation 
@mtkbuild sys = modelingtoolkitize(prob)

prob_jac = ODEProblem(sys, [], (0.0, 1e5), jac=true)
lvjac = eval(generate_jacobian(sys))
plot(sol)