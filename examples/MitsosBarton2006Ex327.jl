using BilevelOptSolver

#include("../src/forrest_solver.jl")
#using .forrest_solver

n₁::Int64 = 5
n₂::Int64 = 5

function F(xy)
    x = @view xy[1:n₁]
    y = @view xy[n₁+1:n₁+n₂]
    x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2 + x[5]^2 + y[1]^2 + y[2]^2 + y[3]^2 + y[4]^2 + y[5]^2
end

function G(xy)
    x = @view xy[1:n₁]
    y = @view xy[n₁+1:n₁+n₂]
    [x[1] + 1; x[2] + 1; x[3] + 1; x[4] + 1; x[5] + 1; 1 - x[1]; 1 - x[2]; 1 - x[3]; 1 - x[4]; 1 - x[5]; x[1] - y[1] * y[2]; exp(x[2]) - y[3] - x[1]; -x[2] * y[1]^2]
end

function f(xy)
    x = @view xy[1:n₁]
    y = @view xy[n₁+1:n₁+n₂]
    y[3] / 10 + y[2]^2 * (x[1] + x[2]) + y[1]^3 + x[3] * x[4] * x[5] * (y[4]^2 + y[5]^2)
end

function g(xy)
    x = @view xy[1:n₁]
    y = @view xy[n₁+1:n₁+n₂]
    [y[1] + 1; y[2] + 1; y[3] + 1; y[4] + 1; y[5] + 1; 1 - y[1]; 1 - y[2]; 1 - y[3]; 1 - y[4]; 1 - y[5]; 3 / 10 - y[1] * y[2]; y[3]^2 - x[1] + 1 / 5; exp(y[3]) - y[4] * y[5] + 1 / 10]
end

xy_init = [1.0; 1; 1; 1; 1; 1; 1; 1; 1; 1]
Ff_optimal = [2.0; -1.1; 2]

bop = construct_bop(n₁, n₂, F, G, f, g; verbosity=0);
sol, is_success, iter_count = solve_bop(bop; x_init=xy_init, verbosity=2, max_iter=200)
if is_success
    @info "success" sol
end

#OP1 = forrest_solver.OptimizationProblem(3, 1:1, F, G, zeros(2), fill(Inf, 2))
#OP2 = forrest_solver.OptimizationProblem(3, 1:2, f, g, zeros(2), fill(Inf, 2))
#bilevel = [OP1; OP2]
#out = forrest_solver.solve(bilevel, [xy_init; zeros(26)])
#sol_forrest = out[1:n]
#@info (sol_forrest)
