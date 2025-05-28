using BilevelOptSolver
using BenchmarkTools
#include("../src/forrest_solver.jl")
#using .forrest_solver

n₁::Int64 = 1
n₂::Int64 = 2

function F(xy)
    x = @view xy[1:n₁]
    y = @view xy[n₁+1:n₁+n₂]
    x[1]
end

function G(xy)
    x = @view xy[1:n₁]
    y = @view xy[n₁+1:n₁+n₂]
    [
        x[1] + 1
        1 - x[1]
    ]
end

function f(xy)
    x = @view xy[1:n₁]
    y = @view xy[n₁+1:n₁+n₂]
    (x[1] - y[1])^2 + (y[2] + 1)^2
end

function g(xy)
    x = @view xy[1:n₁]
    y = @view xy[n₁+1:n₁+n₂]
    [
        y[2] - y[1]^3
        y[2]
    ]
end

xy_init = [1.0; 1; 1]
Ff_optimal = [-1.0; 1; 1]
xy_optimal = [-1.0, -1, 0]

bop = construct_bop(n₁, n₂, F, G, f, g; verbosity=0);
sol, is_success, iter_count = @btime solve_bop(bop; x_init=xy_init, verbosity=0)
if is_success
    @info "success" sol
end

# does not work
#OP1 = forrest_solver.OptimizationProblem(3, 1:1, F, G, zeros(2), fill(Inf, 2))
#OP2 = forrest_solver.OptimizationProblem(3, 1:2, f, g, zeros(2), fill(Inf, 2))
#bilevel = [OP1; OP2]
#out = forrest_solver.solve(bilevel, [xy_init; zeros(26)])
#sol_forrest = out[1:n]
#@info (sol_forrest)
