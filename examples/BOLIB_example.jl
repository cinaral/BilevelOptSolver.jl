# https://biopt.github.io/bolib/
using BilevelOptSolver
include("../src/BOLIB_utils.jl")
using BenchmarkTools
using ProfileView
# Change this to any BOLIB example, e.g.: AllendeStill2013, AnEtal2009, Bard1988Ex2, LamparielloSagratella2017Ex23, Zlobec2001b, AiyoshiShimizu1984Ex2, Colson2002BIPA1

# Quadratic-quadratic Bard1988Ex1, Bard1988Ex2, Bard1988Ex3, Dempe92
# doesn't compile SinhaMaloDeb2014TP9, SinhaMaloDeb2014TP10
b = BOLIB.AnEtal2009()

bop = construct_bop(b.n1, b.n2, b.F, b.G, b.f, b.g, verbosity=0)

elapsed_time = @elapsed begin
    x, status, iter_count = solve_bop(bop; max_iter=200, x_init=b.xy_init, verbosity=5, is_using_HSL=true, is_check_v_agreem=false, tol=1e-6, is_using_PATH_to_init=false, norm_dv_len=10, conv_dv_len=3)
end

is_optimal, is_best, Ff, Ff_star, rating = rate_BOLIB_result(b, bop, x)
print("status: [$status], $iter_count iters ($(round(elapsed_time, sigdigits=5)) s),\t" * rating * ",\t x = $(round.(x, sigdigits=5)) -> Ff = $(round.(Ff, sigdigits=5)) (Ff* = $(round.(Ff_star, sigdigits=5)))\n");

#include("../src/forrest_solver.jl")
#using .forrest_solver
#OP1 = forrest_solver.OptimizationProblem(b.n1 + b.n2, 1:b.n1, b.F, b.G, zeros(bop.m1), fill(Inf, bop.m1))
#OP2 = forrest_solver.OptimizationProblem(b.n1 + b.n2, 1:b.n2, b.f, b.g, zeros(bop.m2), fill(Inf, bop.m2))
#bilevel = [OP1; OP2]
#out = forrest_solver.solve(bilevel)
#x_forrest = out[1:b.n1+b.n2]
#@info x_forrest 
#@info bop.F(x_forrest)