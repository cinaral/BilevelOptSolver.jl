using BilevelOptSolver
# MitsosBarton2006Ex328using BilevelOptSolver

include("../src/forrest_solver.jl")
using .forrest_solver

# Colson2002BIPA1
n₁::Int64 = 1
n₂::Int64 = 1
n::Int64 = n₁ + n₂

function F(xy)
    x = @view xy[1]
    y = @view xy[n₁+1]
    (10 - x[1])^3 + (10 - y[1])^3
end

function G(xy)
    x = @view xy[1]
    y = @view xy[n₁+1]
    [
        5 - x[1]
        x[1] - y[1]
        x[1]
    ]
end

function f(xy)
    x = @view xy[1]
    y = @view xy[n₁+1]
    (x[1] + 2 * y[1] - 15)^4
end

function g(xy)
    x = @view xy[1]
    y = @view xy[n₁+1]
    [
        20 - y[1] - x[1]
        20 - y[1]
        y[1]
    ]
end

Ff_optimal = Float64[250; 0; 1]

x_init = [10; 10]
bop = construct_bop(n₁, n₂, F, G, f, g; verbosity=0);
sol, is_success, iter_count = @btime solve_bop(bop; x_init, verbosity=0)
if is_success
    @info "success" sol
end

OP1 = forrest_solver.OptimizationProblem(2, 1:1, F, G, zeros(3), fill(Inf, 3))
OP2 = forrest_solver.OptimizationProblem(2, 1:1, f, g, zeros(3), fill(Inf, 3))
bilevel = [OP1; OP2]
#out = forrest_solver.solve(bilevel)
out = @btime  forrest_solver.solve(bilevel, [x_init; zeros(33)])
sol_forrest = out[1:n]
@info (sol_forrest)
