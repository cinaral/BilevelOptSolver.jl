# https://biopt.github.io/bolib/
using BilevelOptSolver
include("../src/BOLIB_utils.jl")

# Change this to any BOLIB examples: e.g. AllendeStill2013, AnEtal2009, Bard1988Ex2, LamparielloSagratella2017Ex23, Zlobec2001b, AiyoshiShimizu1984Ex2
"""
Notes 2025-06-23:
AiyoshiShimizu1984Ex2 has Ff* = [5.0, 0.0] but this is wrong?

"""
b = BOLIB.AllendeStill2013()

bop = construct_bop(b.n1, b.n2, b.F, b.G, b.f, b.g, verbosity=0)

elapsed_time = @elapsed begin
    x, is_success, iter_count = solve_bop(bop; max_iter=10, x_init=b.xy_init, verbosity=5, is_using_PATH=false)
end

is_optimal, is_best, Ff, Ff_star, rating = rate_BOLIB_result(b, bop, x)
print("success = $is_success, elapsed: $(round(elapsed_time, sigdigits=5)) s,\t" * rating * "\tx = $(round.(x, sigdigits=5)) -> Ff = $(round.(Ff, sigdigits=5)) (Ff* = $(round.(Ff_star, sigdigits=5)))\n");

#include("../src/forrest_solver.jl")
#using .forrest_solver
#OP1 = forrest_solver.OptimizationProblem(b.n1 + b.n2, 1:b.n1, b.F, b.G, zeros(bop.m₁), fill(Inf, bop.m₁))
#OP2 = forrest_solver.OptimizationProblem(b.n1 + b.n2, 1:b.n2, b.f, b.g, zeros(bop.m₂), fill(Inf, bop.m₂))
#bilevel = [OP1; OP2]
#out = forrest_solver.solve(bilevel)
#x_forrest = out[1:b.n1+b.n2]
#@info x_forrest