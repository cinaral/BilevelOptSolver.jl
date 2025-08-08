# https://biopt.github.io/bolib/
using BilevelOptSolver
using Statistics
using CSV
include("../src/BOLIB_utils.jl")

res = run_all_BOLIB_examples(; verbosity=0, max_iter=50, tol=1e-7, conv_dv_len=3, is_checking_x_agree=true, is_always_hp=false, is_forcing_inactive_inds=false, max_random_restart_count=10, init_solver="PATH", solver="PATH");
CSV.write("BOLIB_results.csv", res.df)

#res_min = run_all_BOLIB_examples(; verbosity=0, max_iter=50, tol=1e-7, conv_dv_len=1, is_always_hp=false, init_solver="PATH", solver="PATH");
#df = hcat(DataFrame("not-minimizing" => fill("", length(res.df[!, 1]))), res.df, DataFrame("minimizing" => fill("", length(res.df[!, 1]))), res_min.df[!, 6:12], makeunique=true)
#CSV.write("BOLIB_results.csv", df)
