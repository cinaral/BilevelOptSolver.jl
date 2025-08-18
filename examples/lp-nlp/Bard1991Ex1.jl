using BilevelOptSolver
""" 
Bard1991Ex1

2025-08-18 works with is_nonstrict_ok=true but not is_nonstrict_ok=false
"""
n1::Int64 = 1
n2::Int64 = 2

function F(xy)
    x = @view xy[1:n1]
    y = @view xy[n1+1:n1+n2]
    x[1] + y[2]
end

function G(xy)
    x = @view xy[1:n1]
    y = @view xy[n1+1:n1+n2]
    [x[1] - 2; 4 - x[1]]
end

function f(xy)
    x = @view xy[1:n1]
    y = @view xy[n1+1:n1+n2]
    2 * y[1] + x[1] * y[2]
end

function g(xy)
    x = @view xy[1:n1]
    y = @view xy[n1+1:n1+n2]
    [y[1] - x[1] + y[2] - 4; y[1]; y[2]]
end

# TODO 2025-08-18: is_checking_x_agree=true leads to bug
xy_optimal = [2.0; 6.0; 0.0]
Ff_optimal = [2.0; 12.0]
x_init = [0.5; -0.5; 1.] 
bop, syms = construct_bop(n1, n2, F, G, f, g; verbosity=0)
is_sol_valid, x, λ, iter_count, status = solve_bop(bop; max_iter=50, x_init, verbosity=5, tol=1e-7, conv_dv_len=3, is_checking_x_agree=false, is_always_hp=false, is_nonstrict_ok=false, max_random_restart_count=100, init_solver="IPOPT", solver="IPOPT")

Ff = [bop.F(x); bop.f(x)]

if is_sol_valid
    print("success ")
else
    print("FAIL ")
end
print("($status), $iter_count iters,\t", "x = $(round.(x, sigdigits=5)) -> Ff = $(round.(Ff, sigdigits=5)) (x* = $(round.(xy_optimal, sigdigits=5)) Ff* = $(round.(Ff_optimal, sigdigits=5)))\n\n")
