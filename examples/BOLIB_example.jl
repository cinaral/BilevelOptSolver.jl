# https://biopt.github.io/bolib/
using BilevelOptSolver
include("../src/BOLIB_utils.jl")

if !haskey(ENV, "BOLIB_PATH")
    error("You need to obtain [BOLIB.jl](https://github.com/cinaral/BOLIB.jl) and set BOLIB_PATH environment variable to run this example.\n")
end
Pkg.develop(path=ENV["BOLIB_PATH"])
import BOLIB

p = BOLIB.Zlobec2001b()

bop = construct_bop(p.n1, p.n2, p.F, p.G, p.f, p.g, verbosity=0)

elapsed_time = @elapsed begin
    x, is_success, iter_count = solve_bop(bop; max_iter=200, x_init=p.xy_init, verbosity=1, is_using_PATH=false)
end

if is_success
    print("success x = $x in $elapsed_time s, Ff = $([bop.F(x); bop.f(x)])\n")
    rate_BOLIB_success(p, bop, x, elapsed_time)
else
    print("fail")
end
