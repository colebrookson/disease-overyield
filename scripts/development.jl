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
    dN₁ = ((ϕ * r₁ᵩ * N₁) * (1 - ϕ) * r₁ * N₁) * (1 - ((N₁ + α₁₂ * N₂) / K₁))
    dN₂ = r₂ * N₂ * (1 - ((N₂ + α₂₁ * N₁) / K₂))
end r₁ r₁ᵩ r₂ α₁₂ α₂₁ K₁ K₂ ϕ
latexify(lv_prop)





































ModelingToolkit.Latexify(lotka_volterra!)

# compute the Jacobian matrix
prob = ODEProblem(lotka_volterra!, [100, 100], (0, 1e5), (1, 1, 0.1, 0.1, 0.2, 0.2))
# modeltoolkit transformation 
@mtkbuild sys = modelingtoolkitize(prob)

prob_jac = ODEProblem(sys, [], (0.0, 1e5), jac=true)
lvjac = eval(generate_jacobian(sys))
plot(sol)