using BilevelOptSolver

# NieEtal2017Ex52
n₁::Int64 = 2
n₂::Int64 = 3

function F(x)
    x₁ = @view x[1:n₁]
    x₂ = @view x[n₁+1:n₁ + n₂]
    x₁[1] * x₂[1] + x₁[2] * x₂[2] + x₁[1] * x₁[2] * x₂[1] * x₂[2] * x₂[3]
end

function G(x)
    x₁ = @view x[1:n₁]
    x₂ = @view x[n₁+1:n₁ + n₂]
    [
        x₁[1] + 1
        x₁[2] + 1
        -x₁[1] + 1
        -x₁[2] + 1
        -x₂[1] * x₂[2] + x₁[1]^2
    ]
end

function f(x)
    x₁ = @view x[1:n₁]
    x₂ = @view x[n₁+1:n₁ + n₂]
    x₁[1] * x₂[1]^2 + x₁[2]^2 * x₂[2] * x₂[3] - x₂[1] * x₂[3]^2
end

function g(x)
    x₁ = @view x[1:n₁]
    x₂ = @view x[n₁+1:n₁ + n₂]
    [
        -1 + x₂[1]^2 + x₂[2]^2 + x₂[3]^2
        -x₂[1]^2 - x₂[2]^2 - x₂[3]^2 + 2
    ]
end

x_init =  [0; 0; 1.1097; 0.3143; −0.8184] # causes a bug
Ff_optimal = [−1; −1; 1.1097; 0.3143; −0.8184]

bop = construct_bop(n₁, n₂, F, G, f, g);
x, is_converged, is_sol_valid, iter_count = solve_bop(bop; x_init, verbosity=2)

Ff = [bop.F(x); bop.f(x)]
print("converged = $is_converged, valid = $is_sol_valid, \tx = $(round.(x, sigdigits=5)) -> Ff = $(round.(Ff, sigdigits=5))\n");

