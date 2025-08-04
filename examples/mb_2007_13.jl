using BilevelOptSolver
using BenchmarkTools
using ProfileView

function mb_2007_13()
    n1::Int64 = 1
    n2::Int64 = 1

    function F(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        x[1] - y[1]
    end

    function G(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        [
            x[1] + 1
            1 - x[1]
        ]
    end

    function f(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        (x[1] * y[1]^2) / 2 - y[1] * x[1]^3
    end

    function g(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        [
            y[1] + 1
            1 - y[1]
        ]
    end

    Ff_optimal = [-1.; 0; 1]
    xy_init = [0.5; 1.]

    (; n1, n2, F, G, f, g, xy_init, Ff_optimal)
end
b = mb_2007_13()
x_optimal = [0.5; 0.5]
#x_init = [0.; 1] # not a solution

bop, syms = construct_bop(b.n1, b.n1, b.F, b.G, b.f, b.g; verbosity=0, np=0)

elapsed_time = @elapsed begin
    is_sol_valid, x, λ, iter_count, status = solve_bop(bop; max_iter=50, x_init=b.xy_init, verbosity=5, tol=1e-7, norm_dv_len=10, conv_dv_len=1, is_checking_min=true, is_checking_x_agree=false, init_solver="IPOPT", solver="PATH")
end

Ff = [bop.F(x); bop.f(x)]
if is_sol_valid
    print("success ")
else
    print("FAIL ")
end
print("($status), $iter_count iters,\t", "x = $(round.(x, sigdigits=5)) -> Ff = $(round.(Ff, sigdigits=5)) (Ff* = $(round.(b.Ff_optimal[1:2], sigdigits=5)))\n");
