
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
function solve_BASBLib_prob(; name="as_2013_01", tol=1e-7, verbosity=5)
    prob = getfield(Main.BASBLib, Symbol(name))()
    bop, syms = construct_bop(prob.n1, prob.n2, prob.F, prob.G, prob.f, prob.g; verbosity=0, np=0)
    x_init = zeros(bop.nx)
    elapsed_time = @elapsed begin
        is_sol_valid, x, λ, iter_count, status = solve_bop(bop; max_iter=50, x_init, verbosity, tol, conv_dv_len=3, do_check_x_agreem=true, do_force_hp_init=false, do_require_strict_min=false, max_rand_restart_ct=50, init_solver="IPOPT", solver="IPOPT")
    end

    Ff = [bop.F(x); bop.f(x)]
    rating = rate_BASBLib_result(name, x, Ff; 1e2 * tol)
    if is_sol_valid
        info_status = "success"
    else
        info_status = "FAIL"
    end
    info = "$info_status: $rating ($status), $iter_count iters ($(round(elapsed_time, sigdigits=5)) s), x = $(round.(x, sigdigits=5)), Ff = $(round.(Ff, sigdigits=5)) (x* = $(round.(prob.xy_optimal, sigdigits=5)), Ff* = $(round.(prob.Ff_optimal, sigdigits=5)))"

    if verbosity > 0
        print(info)
    end

    (; info, Ff, is_sol_valid, x, λ, iter_count, status, elapsed_time, bop, syms)
end

function rate_BASBLib_result(name, x, Ff; tol=1e-7)
    prob = getfield(Main.BASBLib, Symbol(name))()
    if isempty(prob.Ff_optimal)
        return "no sol"
    else
        is_cost_optimal = isapprox(Ff, prob.Ff_optimal; atol=tol)
        is_xy_optimal = isapprox(x, prob.xy_optimal; atol=tol)

        if (is_cost_optimal && is_xy_optimal)
            return "optimal"
        else
            return "NOT optimal"
        end
    end
end

#(; info, Ff, is_sol_valid, x, λ, iter_count, status, elapsed_time, bop, syms) = solve_BASBLib_prob(name="as_2013_01", verbosity=5);
#print(""); # hide output