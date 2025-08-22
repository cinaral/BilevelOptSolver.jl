# https://biopt.github.io/bolib/
using BilevelOptSolver
using Statistics
using CSV
include("../src/BOLIB_utils.jl")

res = run_all_BOLIB_examples(; verbosity=0, max_iter=50, tol=1e-7, conv_dv_len=3, do_check_x_agreem=true, do_force_hp_init=false, max_rand_restart_ct=10, do_require_nonstrict_min=true, init_solver="PATH", solver="PATH");
#CSV.write("BOLIB_results.csv", res.df)

res_min = run_all_BOLIB_examples(; verbosity=0, max_iter=50, tol=1e-7, conv_dv_len=1, do_check_x_agreem=true, do_force_hp_init=false, max_rand_restart_ct=false, init_solver="PATH", solver="PATH");
df = hcat(DataFrame("SOSC" => fill("", length(res.df[!, 1]))), res.df, DataFrame("SONC" => fill("", length(res.df[!, 1]))), res_min.df[!, 6:12], makeunique=true)
CSV.write("BOLIB_results.csv", df)
