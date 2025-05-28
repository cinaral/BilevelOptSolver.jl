using BilevelOptSolver

include("../src/forrest_solver.jl")
using .forrest_solver
using BenchmarkTools

n₁::Int64 = 1
n₂::Int64 = 1
n::Int64 = n₁ + n₂

function F(x)
    x₁ = @view x[1:n₁]
    x₂ = @view x[n₁+1:n]
    x₁[1] + x₂[1]
end

function G(x)
    x₁ = @view x[1:n₁]
    x₂ = @view x[n₁+1:n]
    [
        -x₁[1] - 2 * x₂[1] + 10
        2 * x₁[1] - x₂[1]
        -x₁[1] + 2 * x₂[1]
        x₂[1]
    ]
end

function f(x)
    x₁ = @view x[1:n₁]
    x₂ = @view x[n₁+1:n]
    x₂[1]
end

function g(x)
    x₁ = @view x[1:n₁]
    x₂ = @view x[n₁+1:n]
    [
        x₁[1] + x₂[1] - 3.0,
    ]
end


bop = construct_bop(n₁, n₂, F, G, f, g);
sol, is_success = solve_bop(bop; verbosity=0)
if is_success
    @info "success" sol
end

#OP1 = forrest_solver.OptimizationProblem(2, 1:1, F, G, zeros(4), fill(Inf, 4))
#OP2 = forrest_solver.OptimizationProblem(2, 1:1, f, g, zeros(1), fill(Inf, 1))

#bilevel = [OP1; OP2]
#out = forrest_solver.solve(bilevel)
#sol_forrest = out[1:n]
#@info sol_forrest

