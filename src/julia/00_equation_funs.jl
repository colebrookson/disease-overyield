"""
    dummy_project_function(x, y) → z
Dummy function for illustration purposes.
Performs operation:
```math
z = x + y
```
"""
function dummy_project_function(x, y)
    return x + y
end

function lotka_volterra!(du, u, p, t)
    N1, N2 = du
    r1, r2, K1, α12, α21, K2 = p
    du[1] = dN1 = ((r1 * N1) / K1) * (K1 - N1 - (α12 * N2))
    du[2] = dN2 = ((r2 * N2) / K2) * (K2 - N2 - (α21 * N1))
end
