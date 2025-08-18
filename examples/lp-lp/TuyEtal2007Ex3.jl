using BilevelOptSolver
""" 
TuyEtal2007Ex3

2025-08-18 nothing works
"""

n1::Int64 = 10
n2::Int64 = 6

function F(xy)
    x = @view xy[1:n1]
    y = @view xy[n1+1:n1+n2]
    12 * x[1] - x[2] - 12 * x[3] + 13 * x[4] + 2 * x[6] - 5 * x[8] + 6 * x[9] - 11 * x[10] - 5 * y[1] - 6 * y[2] - 4 * y[3] - 7 * y[4]
end

function G(xy)
    x = @view xy[1:n1]
    y = @view xy[n1+1:n1+n2]
    [2 * x[1] + 3 * x[2] - 14 * x[3] + 2 * x[4] + 9 * x[5] - 2 * x[6] - x[7] + 4 * x[8] - 2 * x[10] + 3 * y[1] - 9 * y[2] + 2 * y[3] + 8 * y[4] - y[5] + 8 * y[6] + 30; 7 * x[2] - x[1] - 13 * x[3] + 15 * x[5] - 2 * x[6] + 8 * x[7] + 4 * x[8] - 4 * x[9] + 7 * x[10] + 6 * y[1] + 2 * y[2] - 6 * y[3] - 2 * y[4] - 8 * y[5] + 4 * y[6] - 134; x[1]; x[2]; x[3]; x[4]; x[5]; x[6]; x[7]; x[8]; x[9]; x[10]]
end

function f(xy)
    x = @view xy[1:n1]
    y = @view xy[n1+1:n1+n2]
    3 * y[1] - 2 * y[2] - 3 * y[3] - 3 * y[4] + y[5] + 6 * y[6]
end

function g(xy)
    x = @view xy[1:n1]
    y = @view xy[n1+1:n1+n2]
    [5 * x[1] - 7 * x[2] - 4 * x[3] + 2 * x[4] - 3 * x[5] + 9 * x[6] - 9 * x[7] + x[8] + 3 * x[9] - 11 * x[10] - 10 * y[1] + 9 * y[2] + 6 * y[3] - 4 * y[4] - 6 * y[5] + 3 * y[6] + 83; 5 * x[2] - 6 * x[1] + 3 * x[3] + 2 * x[4] - 8 * x[5] - 5 * x[6] - 8 * x[7] + 3 * x[8] - 7 * x[9] - 3 * x[10] + 5 * y[1] + 7 * y[2] - y[3] - y[4] + 6 * y[5] - 4 * y[6] + 92; 6 * x[1] + 4 * x[2] - 2 * x[3] + 2 * x[5] - 3 * x[6] + 3 * x[7] - 2 * x[8] - 2 * x[9] - 4 * x[10] - 10 * y[1] - 5 * y[2] - 6 * y[3] + 4 * y[4] - 3 * y[5] + y[6] + 168; 4 * x[4] - 6 * x[2] - 5 * x[1] - 3 * x[5] + 8 * x[6] - x[7] - 2 * x[9] + 3 * x[10] + 4 * y[1] + 3 * y[2] + 4 * y[3] + 4 * y[4] - y[5] - y[6] - 96; 11 * x[2] - 11 * x[1] - 4 * x[3] - 5 * x[4] + 10 * x[5] + 6 * x[6] - 14 * x[7] + 7 * x[8] + 11 * x[9] + 3 * x[10] + 10 * y[1] + 7 * y[2] - 7 * y[3] - 7 * y[4] - 2 * y[5] - 7 * y[6] - 133; 12 * x[2] - 9 * x[1] + 4 * x[3] + 10 * x[4] - 2 * x[5] - 8 * x[6] - 5 * x[7] + 11 * x[8] + 4 * x[9] - x[10] - 2 * y[1] + 5 * y[2] - 10 * y[3] - y[4] - 4 * y[5] - 5 * y[6] + 89; 2 * x[2] - 7 * x[1] + 6 * x[3] + 11 * x[5] - x[6] + 2 * x[7] + 2 * x[8] + x[9] + 2 * x[10] + 5 * y[1] + 5 * y[2] + 6 * y[3] + 5 * y[4] - y[5] + 12 * y[6] - 192; y[1]; y[2]; y[3]; y[4]; y[5]; y[6]]
end

xy_optimal = [0; 8.170692; 10; 0; 7.278940; 3.042311; 0; 10; 0.001982; 9.989153; 3.101280; 10; 10; 10; 0; 9.846133]
Ff_optimal = [-467.4613; -11.6194]

#x_init = [10.0; 10; 10; 10; 10; 10; 10; 10; 10; 10; 10; 10; 10; 10; 10; 10]
x_init = xy_optimal;
bop, syms = construct_bop(n1, n2, F, G, f, g; verbosity=0)
is_sol_valid, x, λ, iter_count, status = solve_bop(bop; max_iter=50, x_init, verbosity=5, tol=1e-7, conv_dv_len=3, is_checking_x_agree=true, is_always_hp=false, is_nonstrict_ok=true, max_random_restart_count=10, init_solver="IPOPT", solver="IPOPT")

Ff = [bop.F(x); bop.f(x)]

if is_sol_valid
    print("success ")
else
    print("FAIL ")
end
print("($status), $iter_count iters,\t", "x = $(round.(x, sigdigits=5)) -> Ff = $(round.(Ff, sigdigits=5)) (x* = $(round.(xy_optimal, sigdigits=5)) Ff* = $(round.(Ff_optimal, sigdigits=5)))\n\n")
