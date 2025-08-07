# https://biopt.github.io/bolib/
using BilevelOptSolver
include("../src/BOLIB_utils.jl")
using BenchmarkTools
using ProfileView
# Change this to any BOLIB example, e.g.: AllendeStill2013, AnEtal2009, Bard1988Ex2, LamparielloSagratella2017Ex23, Zlobec2001b, AiyoshiShimizu1984Ex2, Colson2002BIPA1

# Dimensions 𝑛𝑥, 𝑛𝑦, 𝑛𝐺 or 𝑛𝑔 of examples RobustPortfolioP1, RobustPortfolioP2 can be altered to get problems of different sizes, as necessary

# Zlobec2001b and MitsosBartonEx32 have no optimal solutions
# Quadratic-quadratic Bard1988Ex1, Bard1988Ex2, Bard1988Ex3, Dempe92
# doesn't compile SinhaMaloDeb2014TP9, SinhaMaloDeb2014TP10
# better than optimal AiyoshiShimizu1984Ex2, FalkLiu1995, Mirrlees1999, MitsosBarton2006Ex312, MitsosBarton2006Ex314, MitsosBarton2006Ex315, MitsosBarton2006Ex317, MitsosBarton2006Ex320, PaulaviciusAdjiman2017a, YeZhu2010Ex43, MitsosBarton2006Ex34, MitsosBarton2006Ex35
# FloudasEtal2013 is same as AiyoshiShimizu1984Ex2
"""
better with minimizing: 
LamparielloSagratella2017Ex23 (L-L-N-L) (convex)
GumusFloudas2001Ex1 (N-L-N-L) 
ShimizuEtal1997b (N-L-N-L) 
FrankeEtal2018Ex521 (L-O-L-N) 
Zlobec2001b (no solution but without minimizing we claim there's a local sol) (L-L-O-L-L-N) 

better without minimizing: 
MitsosBarton2006Ex325 (N-N-N-N) 
Outrata1990Ex1c (N-O-N-L) 
AnEtal2009 (N-L-N-L) (convex)
CalamaiVicente1994a (N-O-N-L) 
Dempe1992b (N-O-N-N) (convex)
DempeDutta2012Ex24 (N-O-N-N) 
DempeLohse2011Ex31b (N-O-N-L) 
LamparielloSagratella2017Ex35 (N-L-L-L) 
MitsosBarton2006Ex312 (N-L-N-L) 
NieWangYe2017Ex34 (L-L-N-N) 
FrankeEtal2018Ex513 (L-O-L-N) 
and maybe SinhaMaloDeb2014TP7 (N-N-N-L), TollSettingP2 (N-L-N-L), TollSettingP3 (N-L-N-L)
"""
# interesting ones: AiyoshiShimizu1984Ex2, PaulaviciusAdjiman2017b, Yezza1996Ex41, Mirrlees1999, Outrata1990Ex2e, KleniatiAdjiman2014Ex3
b = BOLIB.KleniatiAdjiman2014Ex4() 

bop, syms = construct_bop(b.n1, b.n2, b.F, b.G, b.f, b.g; verbosity=0, np=0)

x_init = b.xy_init
#x_init = [2.0; -.5]; # gets [1.97; -0.98]
#x_init = [2.0; .5]; # gets [2.1702, 0.4138] 
#x_init = [1.0; 1.]; # gets [1.9912, 0.89471]
x_init = [1; 0.957]; # gets [1.9912, 0.89471]
#x_init = [1.;-.1]
#x_optimal = [1; 0.957]
#x_optimal = [-1.; 1]
elapsed_time = @elapsed begin
    is_sol_valid, x, λ, iter_count, status = solve_bop(bop; max_iter=50, x_init, verbosity=5, tol=1e-7, norm_dv_len=10, conv_dv_len=1, is_checking_min=true, is_checking_x_agree=true, is_always_hp=false, is_forcing_inactive_inds=false, is_require_all_solved=false, init_solver="PATH", solver="PATH")
end

is_optimal, is_best, Ff, Ff_star, rating = rate_BOLIB_result(b, bop, x)
if is_sol_valid
    print("success ")
else
    print("FAIL ")
end
print("($status), $iter_count iters ($(round(elapsed_time, sigdigits=5)) s),\t" * rating * ",\t x = $(round.(x, sigdigits=5)) -> Ff = $(round.(Ff, sigdigits=5)) (Ff* = $(round.(Ff_star, sigdigits=5)))\n");
