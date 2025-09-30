
if !haskey(ENV, "BASBLib_PATH")
    error("You need to obtain [BASBLib.jl](https://github.com/cinaral/BASBLib.jl) (converted from [BASBLibv2.2](https://github.com/basblsolver/BASBLib)) and set BASBLib_PATH environment variable to run this example.\n")
end
Pkg.develop(path=ENV["BASBLib_PATH"])
import BASBLib

using BilevelOptSolver

"""
Usage:
```
julia> include("examples/BASBLib_examples.jl")
julia> (; info, Ff, is_sol_valid, x, λ, iter_count, status, elapsed_time, bop, syms) = solve_BOLIB_prob(name="as_2013_01", verbosity=5);
```
"""
function solve_BASBLib_prob(; name="as_2013_01", x_init=[], tol=1e-7, verbosity=5, init_solver="IPOPT", solver="IPOPT", max_iter=50, conv_dv_len=3, do_force_hp_init=false,  do_check_x_agreem=true, max_rand_restart_ct=10, rating_tol=1e-3)
    prob = getfield(Main.BASBLib, Symbol(name))()
    bop, syms = construct_bop(prob.n1, prob.n2, prob.F, prob.G, prob.f, prob.g; verbosity=0, np=0)

    function solve(x_init=[])
        if isempty(x_init)
            x_init = zeros(bop.nx)
        end
        is_sol_valid, x, λ, iter_count, status, worst_fol_cond, worst_sbop_cond = solve_bop(bop; x_init, verbosity, tol, init_solver, solver, max_iter, conv_dv_len, do_force_hp_init, do_check_x_agreem, max_rand_restart_ct, x_init_min=fill(-10.0, bop.nx), x_init_max=fill(10.0, bop.nx))

        (; is_sol_valid, x, λ, iter_count, status, worst_fol_cond, worst_sbop_cond)
    end

    elapsed_time = @elapsed begin
        is_sol_valid, x, λ, iter_count, status, worst_fol_cond, worst_sbop_cond = solve(x_init)
    end

    Ff = [bop.F(x); bop.f(x)]
    rating = rate_BASBLib_result(name, x, Ff; tol=1e-3)
    if is_sol_valid
        info_status = "success"
    else
        info_status = "FAIL"
    end
    info = "$info_status: $rating ($status $worst_fol_cond $worst_sbop_cond), $iter_count iters ($(round(elapsed_time, sigdigits=5)) s), x = $(round.(x, sigdigits=5)), Ff = $(round.(Ff, sigdigits=5)) (x* = $(round.(prob.xy_optimal, sigdigits=5)), Ff* = $(round.(prob.Ff_optimal, sigdigits=5)))"

    if verbosity > 0
        print(info)
    end

    (; info, Ff, is_sol_valid, x, λ, iter_count, status, worst_fol_cond, worst_sbop_cond, elapsed_time, bop, syms, solve)
end

function rate_BASBLib_result(name, x, Ff; tol=1e-3)
    prob = getfield(Main.BASBLib, Symbol(name))()
    if isempty(prob.Ff_optimal)
        return "no sol"
    else
        is_cost_optimal = isapprox(Ff, prob.Ff_optimal; rtol=tol)
        is_xy_optimal = isapprox(x, prob.xy_optimal; atol=tol)

        if (is_cost_optimal && is_xy_optimal)
            return "optimal"
        else
            return "NOT optimal"
        end
    end
end

#(; info, Ff, is_sol_valid, x, λ, iter_count, status, elapsed_time, bop, syms, solve) = solve_BASBLib_prob(name="as_2013_01", verbosity=0);
#print(""); # hide output

#@btime solve()
#include("../src/forrest_solver.jl")
#using .forrest_solver
#OP1 = forrest_solver.OptimizationProblem(bop.n1 + bop.n2, 1:bop.n1, bop.F, bop.G, zeros(bop.m1), fill(Inf, bop.m1))
#OP2 = forrest_solver.OptimizationProblem(bop.n1 + bop.n2, 1:bop.n2, bop.f, bop.g, zeros(bop.m2), fill(Inf, bop.m2))
#bilevel = [OP1; OP2]
#out_bilevel = @btime forrest_solver.solve(bilevel)
#x_forrest = out_bilevel[1:bop.n1+bop.n2]
#@info x_forrest
#@info bop.F(x_forrest)

#nash = [OP1 OP2]
#out_nash = forrest_solver.solve(nash)