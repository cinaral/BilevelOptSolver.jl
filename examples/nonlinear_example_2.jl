using BilevelOptSolver

include("../src/forrest_solver.jl")
using .forrest_solver
using BenchmarkTools
#using ProfileView

# NieEtal2017Ex52
n₁::Int64 = 2
n₂::Int64 = 3
n::Int64 = n₁ + n₂

function F(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n]
    x[1] * y[1] + x[2] * y[2] + x[1] * x[2] * y[1] * y[2] * y[3]
end

function G(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n]
    [
        x[1] + 1
        x[2] + 1
        -x[1] + 1
        -x[2] + 1
        -y[1] * y[2] + x[1]^2
    ]
end

function f(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n]
    x[1] * y[1]^2 + x[2]^2 * y[2] * y[3] - y[1] * y[3]^2
end

function g(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n]
    [
        -1 + y[1]^2 + y[2]^2 + y[3]^2
        -y[1]^2 - y[2]^2 - y[3]^2 + 2
    ]
end

#x_init =  [0; 0; 1.1097; 0.3143; −0.8184] # causes a bug
global_opt_sol = [−1; −1; 1.1097; 0.3143; −0.8184]

bop = construct_bop(n₁, n₂, F, G, f, g);
#sol = @btime solve_bop(bop; x_init)
sol = @btime solve_bop(bop; x_init=global_opt_sol)

OP1 = forrest_solver.OptimizationProblem(5, 1:2, F, G, zeros(5), fill(Inf, 5))
OP2 = forrest_solver.OptimizationProblem(5, 1:3, f, g, zeros(2), fill(Inf, 2))

bilevel = [OP1; OP2]
out = @btime forrest_solver.solve(bilevel, [x_init; zeros(63)])
#out = @btime forrest_solver.solve(bilevel, [global_opt_sol; zeros(63)])
sol_forrest = out[1:n]

@info sol
@info sol_forrest
