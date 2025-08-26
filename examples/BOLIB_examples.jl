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
function solve_BOLIB_prob(; name="AiyoshiShimizu1984Ex2", tol=1e-7, verbosity=5, rate_fun=rate_BOLIB_result)
    prob = getfield(Main.BOLIB, Symbol(name))()
    bop, syms = construct_bop(prob.n1, prob.n2, prob.F, prob.G, prob.f, prob.g; verbosity=0, np=0)
    x_init = prob.xy_init
    elapsed_time = @elapsed begin
        is_sol_valid, x, λ, iter_count, status = solve_bop(bop; max_iter=50, x_init, verbosity, tol, conv_dv_len=3, do_check_x_agreem=true, do_force_hp_init=false, do_require_nonstrict_min=false, max_rand_restart_ct=50, init_solver="IPOPT", solver="IPOPT")
    end

    Ff = [bop.F(x); bop.f(x)]
    rating = rate_fun(name, x, Ff; tol)
    if is_sol_valid
        info_status = "success"
    else
        info_status = "FAIL"
    end
    info = "$info_status: $rating ($status), $iter_count iters ($(round(elapsed_time, sigdigits=5)) s), x = $(round.(x, sigdigits=5)) -> Ff = $(round.(Ff, sigdigits=5)) (Ff* = $(round.(prob.Ff_optimal, sigdigits=5)))"

    if verbosity > 0
        print(info)
    end

    (; info, Ff, is_sol_valid, x, λ, iter_count, status, elapsed_time, bop, syms)
end

function rate_BOLIB_result(name, x, Ff; tol=1e-2)
    prob = getfield(Main.BOLIB, Symbol(name))()
    is_probably_wrong = false
    is_optimal = false
    is_best = false
    Ff_star = prob.Ff_optimal[1:2]
    rating = ""

    if prob.Ff_optimal[3] == 1 # star means optimal
        if isapprox(Ff, Ff_star; atol=2 * tol)
            is_optimal = true
            rating = "optimal Ff"
        elseif isapprox(Ff[1], Ff_star[1]; atol=2 * tol)
            is_optimal = true
            if Ff[2] < Ff_star[2] - tol
                rating = "optimal F and f is somehow BETTER"
            else
                rating = "optimal F but f is worse"
            end
        elseif Ff[1] < Ff_star[1] - tol || (isapprox(Ff[1], Ff_star[1]; atol=2 * tol) && Ff[2] < Ff_star[2] - tol)
            is_probably_wrong = true
            rating = "better than optimal?!"
        else
            rating = "SUBOPTIMAL Ff"
        end
    elseif prob.Ff_optimal[3] == 2 # star means best
        if isapprox(Ff, Ff_star; atol=2 * tol)
            is_best = true
            rating = "best known Ff"
        elseif isapprox(Ff[1], Ff_star[1]; atol=2 * tol)
            is_best = true
            if Ff[2] < Ff_star[2] - tol
                rating = "best known F and f is somehow BETTER"
            else
                rating = "best known F but f is worse"
            end
        elseif Ff[1] < Ff_star[1] - tol || (isapprox(Ff[1], Ff_star[1]; atol=2 * tol) && Ff[2] < Ff_star[2] - tol)
            is_best = true
            is_probably_wrong = true
            rating = "better than best known Ff"
        else
            rating = "WORSE than best known Ff"
        end
    else
        rating = "no reference solution"
    end
    rating
end

#(; info, Ff, is_sol_valid, x, λ, iter_count, status, elapsed_time, bop, syms) = solve_BOLIB_prob(name="AiyoshiShimizu1984Ex2", verbosity=5);
#print(""); # hide output
