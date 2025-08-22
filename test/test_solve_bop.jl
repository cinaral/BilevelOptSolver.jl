using BilevelOptSolver
using Test
""" 
inspired by [mb_2007_03](https://basblsolver.github.io/BASBLib/LP-QP/mb_2007_03)
convex f
"""

n1::Int64 = 1
n2::Int64 = 1
np::Int64 = 1

function F(xyp)
    x = @view xyp[1:n1]
    y = @view xyp[n1+1:n1+n2]
    p = @view xyp[n1+n2+1:n1+n2+np]
    y[1] - x[1] + p[1]
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

function g(xyp)
    x = @view xyp[1:n1]
    y = @view xyp[n1+1:n1+n2]
    p = @view xyp[n1+n2+1:n1+n2+np]
    [
        y[1]^2 - x[1]
        y[1] + p[1]
        p[1] - y[1]
    ]
end

xy_optimal = [1.0; -1.0]
Ff_optimal = [-1.0; 1.0]

bop, bop_sym = BilevelOptSolver.construct_bop(n1, n2, F, G, f, g; np, verbosity=0)

p::Vector{Float64} = [1.0]
tol::Float64 = 1e-7

# IPOPT:
#x_init = [0.5; -0.5]  # leads to [1;-1] (optimal)
#x_init = [-.5; -.5]  # leads to [1;-1] (optimal)
#x_init = [.5; .5]  # leads to [1;-1] (optimal)
#x_init = [-.5; .5]  # leads to [1;-1] (optimal)
is_sol_valid, x, λ, iter_count, status = BilevelOptSolver.solve_bop(bop; p, max_iter=50, x_init=[0.5; -0.5], verbosity=0, tol, init_solver="IPOPT", solver="IPOPT")
@test isapprox(x, [1; -1])

is_sol_valid, x, λ, iter_count, status = BilevelOptSolver.solve_bop(bop; p, max_iter=50, x_init=[-0.5; -0.5], verbosity=0, tol, init_solver="IPOPT", solver="IPOPT")
@test isapprox(x, [1; -1])

is_sol_valid, x, λ, iter_count, status = BilevelOptSolver.solve_bop(bop; p, max_iter=50, x_init=[0.5; 0.5], verbosity=0, tol, init_solver="IPOPT", solver="IPOPT")
@test isapprox(x, [1; -1])

is_sol_valid, x, λ, iter_count, status = BilevelOptSolver.solve_bop(bop; p, max_iter=50, x_init=[-0.5; 0.5], verbosity=0, tol, init_solver="IPOPT", solver="IPOPT")
@test isapprox(x, [1; -1])

# PATH:
#x_init = [0.5; -0.5]  # leads to [1;-1] (optimal)
#x_init = [-.5; -.5]  # leads to [0;0] (feasible)
#x_init = [.5; .5]  # leads to [0.25;0.5] (feasible)
#x_init = [-.5; .5]  # leads to [0.;0.] (feasible)
is_sol_valid, x, λ, iter_count, status = BilevelOptSolver.solve_bop(bop; p, max_iter=50, x_init=[0.5; -0.5], verbosity=0, tol, init_solver="PATH", solver="PATH")
@test isapprox(x, [1; -1])

is_sol_valid, x, λ, iter_count, status = BilevelOptSolver.solve_bop(bop; p, max_iter=50, x_init=[-0.5; -0.5], verbosity=0, tol, init_solver="PATH", solver="PATH")
@test isapprox(x, [0; 0])

is_sol_valid, x, λ, iter_count, status = BilevelOptSolver.solve_bop(bop; p, max_iter=50, x_init=[0.5; 0.5], verbosity=0, tol, init_solver="PATH", solver="PATH")
@test isapprox(x, [0.25; 0.5])

is_sol_valid, x, λ, iter_count, status = BilevelOptSolver.solve_bop(bop; p, max_iter=50, x_init=[-0.5; 0.5], verbosity=0, tol, init_solver="PATH", solver="PATH")
@test isapprox(x, [0.0; 0.0])

#Ff = [bop.F([x; p]); bop.f([x; p])]

#if is_sol_valid
#    print("success ")
#else
#    print("FAIL ")
#end
#print("($status), $iter_count iters,\t", "x = $(round.(x, sigdigits=5)) -> Ff = $(round.(Ff, sigdigits=5)) (x* = $(round.(xy_optimal, sigdigits=5)) Ff* = $(round.(Ff_optimal, sigdigits=5)))\n\n")
