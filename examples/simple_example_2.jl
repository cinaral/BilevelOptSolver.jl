using BilevelOptSolver

include("../src/forrest_solver.jl")
using .forrest_solver
using BenchmarkTools

# Kleinert (2021)
n₁::Int64 = 1
n₂::Int64 = 1
n::Int64 = n₁ + n₂

function F(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n]
    x[1] + 6y[1]
end

function G(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n]
    [
        x[1] - 5y[1] + 12.5,
    ]
end

function f(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n]
    -y[1]
end

function g(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n]
    [
        2x[1] - y[1]
        -x[1] - y[1] + 6
        -x[1] + 6y[1] + 3
        x[1] + 3y[1] - 3
    ]
end


"""
TODO 2025-05-11:  
1. merge into repo and clear commits DONE
2. fix bop type instability DONE
3. fix feas LP... DONE
4. test on BOLIP
5. make modular
"""

bop = construct_bop(n₁, n₂, F, G, f, g);
#sol = @btime solve_bop(bop)
sol = @btime solve_bop(bop; x_init=[2.4564338234981746; 0.9845259227566776]) # bilevel feasible init
#sol = solve_bop(bop; x_init=[3/7; 6/7]) # optimal init

OP1 = forrest_solver.OptimizationProblem(2, 1:1, F, G, zeros(1), fill(Inf, 1))
OP2 = forrest_solver.OptimizationProblem(2, 1:1, f, g, zeros(4), fill(Inf, 4))

bilevel = [OP1; OP2]
###sol_fmsu = forrest_solver.solve(bilevel) # doesn't work
sol_forrest = @btime forrest_solver.solve(bilevel, [2.4564338234981746; 0.9845259227566776; zeros(37)]) # bilevel feasible init
##sol_fmsu = forrest_solver.solve(bilevel, [3/7; 6/7; zeros(37)]) # optimal init

@info sol
@info sol_forrest[1:n]
