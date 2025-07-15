# https://biopt.github.io/bolib/
using BilevelOptSolver
using Statistics
using CSV
include("../src/BOLIB_utils.jl")

res = run_all_BOLIB_examples(; verbosity=0, max_iter=50, tol=1e-7, is_checking_min=false, conv_dv_len=2, init_solver="IPOPT", solver="IPOPT");
#CSV.write("BOLIB_results.csv", res.df)

res_min = run_all_BOLIB_examples(; verbosity=0, max_iter=50, tol=1e-7, is_checking_min=true, conv_dv_len=2, init_solver="IPOPT", solver="IPOPT");
df = hcat(DataFrame("PATH & not-minimizing" => fill("", length(res.df[!, 1]))), res.df, DataFrame("IPOPT & minimizing" => fill("", length(res.df[!, 1]))), res_min.df[!, 6:12], makeunique=true)
CSV.write("BOLIB_results.csv", df)
