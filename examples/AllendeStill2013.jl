using BilevelOptSolver
# MitsosBarton2006Ex328using BilevelOptSolver

include("../src/forrest_solver.jl")
using .forrest_solver
using BenchmarkTools
using Random
#using ProfileView

# AllendeStill2013
n₁::Int64 = 2
n₂::Int64 = 2
n::Int64 = n₁ + n₂

function F(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n]
    x[1]^2 - 2 * x[1] + x[2]^2 - 2*x[2] + y[1]^2 + y[2]^2
end

function G(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n]
    [
		x[1]
		x[2]
		y[1]
		y[2]
		-x[1] + 2
    ]
end

function f(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n]
	y[1]^2 - 2*x[1]*y[1] + y[2]^2 - 2*x[2]*y[2]
end

function g(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n]
    [
		-(y[1] - 1)^2 + 0.25
		-(y[2] - 1)^2 + 0.25
    ]
end

x_optimal = [0.5 0.5 0.5 0.5]

bop = construct_bop(n₁, n₂, F, G, f, g; verbosity=2);
sol, _ = solve_bop(bop; verbosity=2)
@info (sol)
#sol = solve_bop(bop; x_init)

OP1 = forrest_solver.OptimizationProblem(4, 1:2, F, G, zeros(5), fill(Inf, 5))
OP2 = forrest_solver.OptimizationProblem(4, 1:2, f, g, zeros(2), fill(Inf, 2))
bilevel = [OP1; OP2]
out = forrest_solver.solve(bilevel, [x_init; zeros(145)])
sol_forrest = out[1:n]
@info (sol_forrest)
