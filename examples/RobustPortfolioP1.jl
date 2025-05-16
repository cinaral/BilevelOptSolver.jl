using BilevelOptSolver

include("../src/forrest_solver.jl")
using .forrest_solver
#using BenchmarkTools
#using ProfileView

n₁::Int64 = 11
n₂::Int64 = 10
n::Int64 = n₁ + n₂

function F(xy)
    x = @view xy[1:n₁]
    y = @view xy[n₁+1:n₁+n₂]
    -x[11]
end

function G(xy)
    x = @view xy[1:n₁]
    y = @view xy[n₁+1:n₁+n₂]
    [
        x[1] * y[1] - x[11] + x[2] * y[2] + x[3] * y[3] + x[4] * y[4] + x[5] * y[5] + x[6] * y[6] + x[7] * y[7] + x[8] * y[8] + x[9] * y[9] + x[10] * y[10]
        x[1]
        x[2]
        x[3]
        x[4]
        x[5]
        x[6]
        x[7]
        x[8]
        x[9]
        x[10]
        1 - x[2] - x[3] - x[4] - x[5] - x[6] - x[7] - x[8] - x[9] - x[10] - x[1]
        x[1] + x[2] + x[3] + x[4] + x[5] + x[6] + x[7] + x[8] + x[9] + x[10] - 1
    ]
end

function f(xy)
    x = @view xy[1:n₁]
    y = @view xy[n₁+1:n₁+n₂]
    x[1] * y[1] - x[11] + x[2] * y[2] + x[3] * y[3] + x[4] * y[4] + x[5] * y[5] + x[6] * y[6] + x[7] * y[7] + x[8] * y[8] + x[9] * y[9] + x[10] * y[10]
end

function g(xy)
    x = @view xy[1:n₁]
    y = @view xy[n₁+1:n₁+n₂]
    [
        9 / 4 - (9000 * abs(y[2] - 29 / 25)^2) / 11 - (3600 * abs(y[5] - 47 / 40)^2) / 11 - (3000 * abs(y[6] - 59 / 50)^2) / 11 - (4500 * abs(y[4] - 117 / 100)^2) / 11 - (2250 * abs(y[8] - 119 / 100)^2) / 11 - (18000 * abs(y[1] - 231 / 200)^2) / 11 - (6000 * abs(y[3] - 233 / 200)^2) / 11 - (18000 * abs(y[7] - 237 / 200)^2) / 77 - (2000 * abs(y[9] - 239 / 200)^2) / 11 - (1800 * abs(y[10] - 6 / 5)^2) / 11
        y[1]
        y[2]
        y[3]
        y[4]
        y[5]
        y[6]
        y[7]
        y[8]
        y[9]
        y[10]
    ]
end

xy_init = [1.0; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1]
Ff_best_known = [1.15; 0; 2]


bop = construct_bop(n₁, n₂, F, G, f, g; verbosity=0);
sol, is_success, iter_count = solve_bop(bop; x_init=xy_init, verbosity=2, max_iter=200)
if is_success
    @info "success" sol
end

#OP1 = forrest_solver.OptimizationProblem(21, 1:11, F, G, zeros(13), fill(Inf, 13))
#OP2 = forrest_solver.OptimizationProblem(21, 1:10, f, g, zeros(11), fill(Inf, 11))
#bilevel = [OP1; OP2]
#out = @btime forrest_solver.solve(bilevel, [xy_init; zeros(144)])
#sol_forrest = out[1:n]
#@info (sol_forrest)
