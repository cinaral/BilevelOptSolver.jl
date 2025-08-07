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

#include("../tests/lp-qp/as_1984_01.jl")
#p = as_1984_01()
##x_init = [3.5449013142516694, 2.4259541458606315, 9.571562944376396, 7.4268034008199955]
##x_init = [0.0; 30; -10; 10]
##x_init = [23.518569764623052, 2.5581337650600555, 1.6513845230341495, 12.191872247621625]
#x_init = [33.625734921878994, 6.882919488425432, 14.884499926602707, 14.05110377850071]
#bop, syms = construct_bop(p.n1, p.n2, p.F, p.G, p.f, p.g; verbosity=0, np=0)
##x_init = [0.0; 0; -10.; -10.] .+ rand(MersenneTwister(123), p.n1 + p.n2) .* [50.0; 50.0; 30.0; 30.0]
#is_sol_valid, x, λ, iter_count, status = solve_bop(bop; max_iter=50, x_init, verbosity=5, tol=1e-7, norm_dv_len=10, conv_dv_len=1, is_checking_min=true, is_checking_x_agree=true, is_forcing_inactive_inds=false, is_always_hp=false, init_solver="IPOPT", solver="IPOPT")

#Ff = [bop.F(x); bop.f(x)]

#if is_sol_valid
#    print("success ")
#else
#    print("FAIL ")
#end
#print("($status), $iter_count iters,\t", "x = $(round.(x, sigdigits=5)) -> Ff = $(round.(Ff, sigdigits=5)) (x* = $(round.(p.xy_optimal, sigdigits=5)) -> Ff* = $(round.(p.Ff_optimal, sigdigits=5)))\n\n")

#include("../tests/qp-qp/as_1981_01.jl")
#p = as_1981_01()
##x_init=[9.103373216185842, 2.3203358259985274, 10.425466604246346, 7.968079005512396, 7.632104096637655, 16.333910141808303, 36.19995643558407, 12.828408932354609]
#x_init = [2.891288336896034, 3.153294421248287, 7.145553084160766, 5.657353341126381, 17.732454127949893, 10.237586457164474, 20.131480429252385, 7.744882067431531]
#bop, syms = construct_bop(p.n1, p.n2, p.F, p.G, p.f, p.g; verbosity=0, np=0)
##x_init = randn(MersenneTwister(123), p.n1 + p.n2)
#is_sol_valid, x, λ, iter_count, status = solve_bop(bop; max_iter=50, x_init, verbosity=5, tol=1e-7, norm_dv_len=10, conv_dv_len=1, is_checking_min=true, is_checking_x_agree=true, is_forcing_inactive_inds=false, is_always_hp=false, init_solver="IPOPT", solver="IPOPT")

#Ff = [bop.F(x); bop.f(x)]

#if is_sol_valid
#    print("success ")
#else
#    print("FAIL ")
#end
#print("($status), $iter_count iters,\t", "x = $(round.(x, sigdigits=5)) -> Ff = $(round.(Ff, sigdigits=5)) (x* = $(round.(p.xy_optimal, sigdigits=5)) -> Ff* = $(round.(p.Ff_optimal, sigdigits=5)))\n\n")