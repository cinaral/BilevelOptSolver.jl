using BilevelOptSolver

#include("../src/forrest_solver.jl")
#using .forrest_solver
#using BenchmarkTools
#using ProfileView

# AiyoshiShimizu1984Ex2
n₁::Int64 = 2
n₂::Int64 = 2
n::Int64 = n₁ + n₂

function F(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n]
    2x[1] + 2x[2] - 3y[1] - 3y[2] - 60
end

function G(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n]
    [
        -x[1] - x[2] - y[1] + 2y[2] + 40
        -x[1] + 50
        -x[2] + 50
        x[1]
        x[2]
    ]
end

function f(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n]
    (y[1] - x[1] + 20)^2 + (y[2] - x[2] + 20)^2
end

function g(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n]
    [
        -2y[1] + x[1] - 10
        -2y[2] + x[2] - 10
        y[1] + 10
        y[2] + 10
        -y[1] + 20
        -y[2] + 20
    ]
end

global_opt_sol = [25.; 30; 5; 10]
local_opt_sol = [0.; 0; -10; -10]

bop = construct_bop(n₁, n₂, F, G, f, g);
sol, is_success, iter_count = solve_bop(bop; verbosity=2)
#sol = solve_bop(bop; x_init=global_opt_sol)
if is_success
    @info "success" sol
end

#OP1 = forrest_solver.OptimizationProblem(4, 1:2, F, G, zeros(5), fill(Inf, 5))
#OP2 = forrest_solver.OptimizationProblem(4, 1:2, f, g, zeros(6), fill(Inf, 6))
#bilevel = [OP1; OP2]
#out = forrest_solver.solve(bilevel)
##out = @btime forrest_solver.solve(bilevel, [global_opt_sol; zeros(64)])
#sol_forrest = out[1:n]
#@info sol_forrest