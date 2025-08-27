using BilevelOptSolver
""" 
inspired by [mb_2007_03](https://basblsolver.github.io/BASBLib/LP-QP/mb_2007_03)
convex f
"""

n1::Int64 = 1
n2::Int64 = 1

function F(xy)
    x = @view xy[1:n1]
    y = @view xy[n1+1:n1+n2]
    y[1] - x[1] + 1.0
end

function G(xy)
    x = @view xy[1:n1]
    y = @view xy[n1+1:n1+n2]
    [
        x[1] + 1.0
        1.0 - x[1]
        y[1] + 1.0
        1.0 - y[1]
    ]
end

function f(xy)
    x = @view xy[1:n1]
    y = @view xy[n1+1:n1+n2]
    y[1]^2
end

function g(xy)
    x = @view xy[1:n1]
    y = @view xy[n1+1:n1+n2]
    [
        y[1]^2 - x[1]
        y[1] + 1.0
        1.0 - y[1]
    ]
end

xy_optimal = [1.0; -1.0]
Ff_optimal = [-1.0; 1.0]

x_init = [0.5; -0.5]  # leads to [1;-1] (optimal)
#x_init = [-.5; -.5]  # leads to [0;0] (feasible)
#x_init = [.5; .5]  # leads to [0.25;0.5] (feasible)
#x_init = [-.5; .5]  # leads to [0.;0.] (feasible)
bop, syms = construct_bop(n1, n2, F, G, f, g; verbosity=0)
param = [1.0; 1.0]
is_sol_valid, x, λ, iter_count, status = solve_bop(bop; max_iter=50, x_init, verbosity=5, tol=1e-7, conv_dv_len=3, do_check_x_agreem=true, do_force_hp_init=false, do_require_strict_min=false, max_rand_restart_ct=50, init_solver="IPOPT", solver="IPOPT")

Ff = [bop.F([x; param]); bop.f([x; param])]

if is_sol_valid
    print("success ")
else
    print("FAIL ")
end
print("($status), $iter_count iters,\t", "x = $(round.(x, sigdigits=5)) -> Ff = $(round.(Ff, sigdigits=5)) (x* = $(round.(xy_optimal, sigdigits=5)) Ff* = $(round.(Ff_optimal, sigdigits=5)))\n")
