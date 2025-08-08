using BilevelOptSolver
using BenchmarkTools
using ProfileView

n1::Int64 = 1
n2::Int64 = 1

function F(x)
    x₁ = @view x[1:n1]
    x₂ = @view x[n1+1:n1+n2]
    #-x₁[1]^2
    (x₁[1] - 1)^2 + (x₂[1] - 1)^2
end

function G(x)
    x₁ = @view x[1:n1]
    x₂ = @view x[n1+1:n1+n2]
    [
    ]
end

function f(x)
    x₁ = @view x[1:n1]
    x₂ = @view x[n1+1:n1+n2]
    -x₂[1]
end

function g(x)
    x₁ = @view x[1:n1]
    x₂ = @view x[n1+1:n1+n2]
    [
        -x₁[1] - x₂[1] + 1
        x₁[1] - x₂[1] + 1
    ]
end
Ff_optimal = [0.5; -0.5; 1]
x_optimal = [0.5; 0.5]

bop, syms = construct_bop(n1, n2, F, G, f, g; verbosity=0, np=0)

#x_init = [-0.1; 1] # leads to [0.; 1] # not a solution
x_init = [0.1; 1] # leads to [0.5; .5] # optimal solution

is_sol_valid, x, λ, iter_count, status = solve_bop(bop; max_iter=50, x_init, verbosity=5, tol=1e-7, conv_dv_len=1, is_checking_x_agree=false, init_solver="PATH", solver="PATH")

Ff = [bop.F(x); bop.f(x)]
if is_sol_valid
    print("success ")
else
    print("FAIL ")
end
print("($status), $iter_count iters,\t", "x = $(round.(x, sigdigits=5)) -> Ff = $(round.(Ff, sigdigits=5)) (Ff* = $(round.(b.Ff_optimal[1:2], sigdigits=5)))\n");
