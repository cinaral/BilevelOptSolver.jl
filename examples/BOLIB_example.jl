# https://biopt.github.io/bolib/
using BilevelOptSolver
include("../src/BOLIB_utils.jl")
using BenchmarkTools
using ProfileView
# Change this to any BOLIB example, e.g.: AllendeStill2013, AnEtal2009, Bard1988Ex2, LamparielloSagratella2017Ex23, Zlobec2001b, AiyoshiShimizu1984Ex2, Colson2002BIPA1

# Dimensions 𝑛𝑥, 𝑛𝑦, 𝑛𝐺 or 𝑛𝑔 of examples RobustPortfolioP1, RobustPortfolioP2 can be altered to get problems of different sizes, as necessary
"""
Interesting ones (F fails, B better than optimal)
---
AnEtal2009 F

DempeDutta2012Ex24 F
DempeDutta2012Ex31 F
FalkLiu1995 B
GumusFloudas2001Ex1 F
HenrionSurowiec2011 F
KleniatiAdjiman2014Ex3 B
LamparielloSagratella2017Ex32 F
LuDebSinha2016f F
MacalHurter1997 F
MitsosBarton2006Ex38 F
MitsosBarton2006Ex39 F
MitsosBarton2006Ex314 B
MitsosBarton2006Ex315 B
MitsosBarton2006Ex318 F 
MitsosBarton2006Ex320 B
MitsosBarton2006Ex321 F
MitsosBarton2006Ex324 F
MitsosBarton2006Ex327 B
NieWangYe2017Ex34 F
Outrata1990Ex1e B
Outrata1990Ex2d B
ShimizuEtal1997b F
WanWangLv2011 B
YeZhu2010Ex42 F
Zlobec2001b F (no solution exists)
DesignCentringP2 F
NetworkDesignP1 F
NetworkDesignP2 F
RobustPortfolioP2 F
TuyEtal2007Ex3 F
MitsosBarton2006Ex32 F (no solution exists)
MitsosBarton2006Ex35 B
MitsosBarton2006Ex36 F
"""
b = BOLIB.TuyEtal2007Ex3()

bop, syms = construct_bop(b.n1, b.n2, b.F, b.G, b.f, b.g; verbosity=0, np=0)
x_init = b.xy_init
elapsed_time = @elapsed begin
    is_sol_valid, x, λ, iter_count, status = solve_bop(bop; max_iter=50, x_init, verbosity=5, tol=1e-7, conv_dv_len=3, is_checking_x_agree=true, is_always_hp=false, is_forcing_inactive_inds=false, is_nonstrict_ok=false, max_random_restart_count=50, init_solver="PATH", solver="PATH")
end

is_optimal, is_best, Ff, Ff_star, rating = rate_BOLIB_result(b, bop, x)
if is_sol_valid
    print("success ")
else
    print("FAIL ")
end
print("($status), $iter_count iters ($(round(elapsed_time, sigdigits=5)) s),\t" * rating * ",\t x = $(round.(x, sigdigits=5)) -> Ff = $(round.(Ff, sigdigits=5)) (Ff* = $(round.(Ff_star, sigdigits=5)))\n");
