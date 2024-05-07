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
    r1, r2, α11, α12, α21, α22 = p
    du[1] = dN1 = N1 * (r1 - (α11 * N1) - (α12 * N2))
    du[2] = dN2 = N2 * (r2 - (α11 * N1) - (α22 * N2))
end

