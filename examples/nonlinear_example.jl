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
x_optimal = [−1; −1; 1.1097; 0.3143; −0.8184]

bop, _ = construct_bop(n₁, n₂, F, G, f, g);

is_sol_valid, x, λ, iter_count, status = solve_bop(bop; x_init, solver="PATH", verbosity=5)

Ff_optimal = [bop.F(x_optimal); bop.f(x_optimal)]
Ff = [bop.F(x); bop.f(x)]
if is_sol_valid
    print("success ")
else
    print("FAIL ")
end
print("($status), $iter_count iters,\t", "x = $(round.(x, sigdigits=5)) -> Ff = $(round.(Ff, sigdigits=5)) (Ff* = $(round.(Ff_optimal, sigdigits=5)))\n");
