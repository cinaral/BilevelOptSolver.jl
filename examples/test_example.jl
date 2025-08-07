using BilevelOptSolver
using Random
include("../tests/lp-lp/bf_1982_02.jl")
p = bf_1982_02()
bop, syms = construct_bop(p.n1, p.n2, p.F, p.G, p.f, p.g; verbosity=0, np=0)
x_init = randn(MersenneTwister(456), p.n1 + p.n2)
is_sol_valid, x, λ, iter_count, status = solve_bop(bop; max_iter=50, x_init, verbosity=5, tol=1e-5, norm_dv_len=10, conv_dv_len=1, is_checking_min=true, is_checking_x_agree=true, is_always_hp=false, init_solver="IPOPT", solver="IPOPT")

Ff = [bop.F(x); bop.f(x)]

if is_sol_valid
    print("success ")
else
    print("FAIL ")
end
print("($status), $iter_count iters,\t", "x = $(round.(x, sigdigits=5)) -> Ff = $(round.(Ff, sigdigits=5)) (x* = $(round.(p.xy_optimal, sigdigits=5)) -> Ff* = $(round.(p.Ff_optimal, sigdigits=5)))\n\n")


#include("../tests/lp-lp/ct_1982_01.jl")
#p = ct_1982_01()
#bop, syms = construct_bop(p.n1, p.n2, p.F, p.G, p.f, p.g; verbosity=0, np=0)
#x_init = randn(MersenneTwister(456), p.n1 + p.n2)
#is_sol_valid, x, λ, iter_count, status = solve_bop(bop; max_iter=50, x_init, verbosity=5, tol=1e-7, norm_dv_len=10, conv_dv_len=1, is_checking_min=true, is_checking_x_agree=true, is_forcing_inactive_inds=false, is_always_hp=false, init_solver="IPOPT", solver="IPOPT")

#Ff = [bop.F(x); bop.f(x)]

#if is_sol_valid
#    print("success ")
#else
#    print("FAIL ")
#end
#print("($status), $iter_count iters,\t", "x = $(round.(x, sigdigits=5)) -> Ff = $(round.(Ff, sigdigits=5)) (x* = $(round.(p.xy_optimal, sigdigits=5)) -> Ff* = $(round.(p.Ff_optimal, sigdigits=5)))\n\n")
