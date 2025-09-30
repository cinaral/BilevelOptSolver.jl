if !haskey(ENV, "BOLIB_PATH")
    error("You need to obtain [BOLIB.jl](https://github.com/cinaral/BOLIB.jl) (converted from [BOLIBver2](https://biopt.github.io/bolib/)) and set BOLIB_PATH environment variable to run this example.\n")
end
Pkg.develop(path=ENV["BOLIB_PATH"])
import BOLIB
using BilevelOptSolver

"""
Usage:
```
julia> include("examples/BOLIB_examples.jl")
julia> (; info, Ff, is_sol_valid, x, λ, iter_count, status, elapsed_time, bop, syms) = solve_BOLIB_prob(name="AiyoshiShimizu1984Ex2", verbosity=5);
```
"""
function solve_BOLIB_prob(; name="AiyoshiShimizu1984Ex2", x_init = [], tol=1e-7, verbosity=5, init_solver="IPOPT", solver="IPOPT", max_iter=50, conv_dv_len=3, do_force_hp_init=false, do_check_x_agreem=true, max_rand_restart_ct=10, rating_tol=1e-3)
    prob = getfield(Main.BOLIB, Symbol(name))()
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
    rating = rate_BOLIB_result(name, x, Ff; tol=rating_tol)
    if is_sol_valid
        info_status = "success"
    else
        info_status = "FAIL"
    end
    info = "$info_status: $rating ($status $worst_fol_cond $worst_sbop_cond), $iter_count iters ($(round(elapsed_time, sigdigits=5)) s), x = $(round.(x, sigdigits=5)), Ff = $(round.(Ff, sigdigits=5)) (Ff* = $(round.(prob.Ff_optimal, sigdigits=5)))"

    if verbosity > 0
        print(info)
    end

    (; info, Ff, is_sol_valid, x, λ, iter_count, status, elapsed_time, bop, syms)
end

function rate_BOLIB_result(name, x, Ff; tol=1e-3)
    prob = getfield(Main.BOLIB, Symbol(name))()

    if prob.Ff_optimal[3] == 0 # star means optimal
        rating = "no reference"
    else
        is_cost_optimal = isapprox(Ff, prob.Ff_optimal[1:2]; rtol=tol)
        if is_cost_optimal
            return "optimal"
        else
            return "NOT optimal"
        end
    end
    rating
end

#(; info, Ff, is_sol_valid, x, λ, iter_count, status, elapsed_time, bop, syms) = solve_BOLIB_prob(name="AiyoshiShimizu1984Ex2", verbosity=5);
#print(""); # hide output
