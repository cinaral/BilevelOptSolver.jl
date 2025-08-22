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

x_init = [0.5; -0.5]  # leads to [1;-1] (optimal)
#x_init = [-.5; -.5]  # leads to [0;0] (feasible)
#x_init = [.5; .5]  # leads to [0.25;0.5] (feasible)
#x_init = [-.5; .5]  # leads to [0.;0.] (feasible)
v = [x_init; zeros(bop.fol_nlp.m)]

BilevelOptSolver.initialize_z!(v, bop; p=1.0, verbosity=0, init_solver="PATH", tol=1e-6, max_iter=100)
@test isapprox(v, [0.5; -0.7071067811873449; 1.0; 0; 0])

BilevelOptSolver.initialize_z!(v, bop; p=1.0, verbosity=0, init_solver="PATH", tol=1e-6, max_iter=100)
@test isapprox(v, [0.5; -0.7071067811873449; 1.0; 0; 0])

#BilevelOptSolver.initialize_z!(v, bop; p=1.0, verbosity=0, init_solver="IPOPT", tol=1e-6, max_iter=100, is_always_hp=true)

#bop, bop_sym = BilevelOptSolver.construct_bop_2(n1, n2, F, G, f, g; verbosity=0)