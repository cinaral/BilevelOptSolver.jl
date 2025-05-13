# MitsosBarton2006Ex328using BilevelOptSolver

include("../src/forrest_solver.jl")
using .forrest_solver
using BenchmarkTools
#using ProfileView

# NieEtal2017Ex52
n₁::Int64 = 5
n₂::Int64 = 5
n::Int64 = n₁ + n₂

function F(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n]
    -(x[1]^2 + y[1]^2) - (x[2]^2 + y[2]^2) - (x[3]^2 + y[3]^2) - (x[4]^2 + y[4]^2) - (x[5]^2 + y[5]^2)
end

function G(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n]
    [
        x[1] + 1
		x[2] + 1
		x[3] + 1
		x[4] + 1
        x[5] + 1
		-x[1] + 1
		-x[2] + 1
		-x[3] + 1
		-x[4] + 1
        -x[5] + 1
        -y[1] * y[2] + x[1]
		-x[2] * y[1]^2
		-x[1] + exp(x[2]) - y[3]
    ]
end

function f(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n]
   y[1]^3 + y[2]^2 * x[1] + y[2]^2 * x[2] + 0.1 * y[3] + (y[4]^2 + y[5]^2) * x[3] * x[4] * x[5]
end

function g(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n]
    [
		y[1] + 1
		y[2] + 1
		y[3] + 1
		y[4] + 1
        y[5] + 1
		-y[1] + 1
		-y[2] + 1
		-y[3] + 1
		-y[4] + 1
        -y[5] + 1
		-y[1] * y[2] + 0.3
		-x[1] + y[3]^2 + 0.2
		exp(y[3]) - y[4] * y[5] + 0.1
    ]
end

#x_init =  [0; 0; 1.1097; 0.3143; −0.8184] # causes a bug
best_x = [1.; −1; -1; −1; -1; -1; 1; −1; −1; 1]
x_init = best_x .+ 1 .- 2 .* rand(MersenneTwister(12345), 10)

bop = construct_bop(n₁, n₂, F, G, f, g);
#sol = solve_bop(bop; x_init=best_x)
#sol = solve_bop(bop; max_iter=3000)
sol = solve_bop(bop; x_init)

OP1 = forrest_solver.OptimizationProblem(10, 1:5, F, G, zeros(13), fill(Inf, 13))
OP2 = forrest_solver.OptimizationProblem(10, 1:5, f, g, zeros(13), fill(Inf, 13))

bilevel = [OP1; OP2]
out = forrest_solver.solve(bilevel, [x_init; zeros(145)])
sol_forrest = out[1:n]

@info sol
@info sol_forrest
