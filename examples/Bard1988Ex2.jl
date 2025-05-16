using BilevelOptSolver
# MitsosBarton2006Ex328using BilevelOptSolver

include("../src/forrest_solver.jl")
using .forrest_solver

# Bard1988Ex2
n₁::Int64 = 4
n₂::Int64 = 4
n::Int64 = n₁ + n₂

function F(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n]
    -(200 - y[1] - y[3]) * (y[1] + y[3]) - (160 - y[2] - y[4]) * (y[2] + y[4])
end

function G(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n]
    [
        -x[1] - x[2] - x[3] - x[4] + 40
        10 - x[1]
        5 - x[2]
        15 - x[3]
        20 - x[4]
        x[1]
        x[2]
        x[3]
        x[4]
    ]
end

function f(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n]
    (y[1] - 4)^2 + (y[2] - 13)^2 + (y[3] - 35)^2 + (y[4] - 2)^2
end

function g(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n]
    [
        x[1] - (2 * y[1]) / 5 - (7 * y[2]) / 10
        x[2] - (3 * y[1]) / 5 - (3 * y[2]) / 10
        x[3] - (2 * y[3]) / 5 - (7 * y[4]) / 10
        x[4] - (3 * y[3]) / 5 - (3 * y[4]) / 10
        20 - y[1]
        20 - y[2]
        40 - y[3]
        40 - y[4]
        y[1]
        y[2]
        y[3]
        y[4]
    ]
end

x_init = [5.0; 5; 15; 15; 0; 0; 0; 0]
bop = construct_bop(n₁, n₂, F, G, f, g; verbosity=0);
sol, is_success, iter_count = solve_bop(bop; x_init, verbosity=1, max_iter=5000)
if is_success
    @info "success" sol
end

OP1 = forrest_solver.OptimizationProblem(8, 1:4, F, G, zeros(9), fill(Inf, 9))
OP2 = forrest_solver.OptimizationProblem(8, 1:4, f, g, zeros(12), fill(Inf, 12))
bilevel = [OP1; OP2]
out = forrest_solver.solve(bilevel)
out = forrest_solver.solve(bilevel, [x_init; zeros(126)])
sol_forrest = out[1:n]
@info (sol_forrest)
