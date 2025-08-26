
if !haskey(ENV, "BASBLib_PATH")
    error("You need to obtain [BASBLib.jl](https://github.com/cinaral/BASBLib.jl) (converted from [BASBLibv2.2](https://github.com/basblsolver/BASBLib)) and set BASBLib_PATH environment variable to run this example.\n")
end
Pkg.develop(path=ENV["BASBLib_PATH"])
import BASBLib

using BilevelOptSolver

function solve_BASBLib_prob(; name="as_2013_01")
    b = getfield(Main.BASBLib, Symbol(name))()

    bop, syms = construct_bop(b.n1, b.n2, b.F, b.G, b.f, b.g; verbosity=0, np=0)
    x_init = zeros(bop.nx)
    elapsed_time = @elapsed begin
        is_sol_valid, x, λ, iter_count, status = solve_bop(bop; max_iter=50, x_init, verbosity=5, tol=1e-7, conv_dv_len=3, do_check_x_agreem=true, do_force_hp_init=false, do_require_nonstrict_min=false, max_rand_restart_ct=50, init_solver="IPOPT", solver="IPOPT")
    end

    Ff = [bop.F(x); bop.f(x)]

    #is_optimal, is_best, Ff, Ff_star, rating = rate_BOLIB_result(b, bop, x)
    if is_sol_valid
        print("success ")
    else
        print("FAIL ")
    end
    print("($status), $iter_count iters ($(round(elapsed_time, sigdigits=5)) s),\t", "x = $(round.(x, sigdigits=5)) -> Ff = $(round.(Ff, sigdigits=5)) (x* = $(round.(b.xy_optimal, sigdigits=5)) Ff* = $(round.(b.Ff_optimal, sigdigits=5)))\n\n")
end

solve_BASBLib_prob(name="as_2013_01")