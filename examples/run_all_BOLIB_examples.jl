# https://biopt.github.io/bolib/
using BilevelOptSolver
using Statistics
using CSV
include("../src/BOLIB_utils.jl")

#@info "is_checking_min FALSE"
res = run_all_BOLIB_examples(; verbosity=0, max_iter=100, tol=1e-7, is_checking_min=true, conv_dv_len=1, init_solver="IPOPT", solver="IPOPT");
CSV.write("BOLIB_results.csv", res.df)

#@info "is_checking_min TRUE"
#res_no_min = run_all_BOLIB_examples(;verbosity=0, max_iter=100, tol=1e-6, is_checking_min=true, conv_dv_len=1, init_solver="IPOPT", solver="IPOPT");

#df = hcat(res.df, DataFrame("non-minimizing" => fill("", length(res.df[!, 1]))), res_no_min.df[!, 6:12], makeunique=true)
#CSV.write("BOLIB_results.csv", df)

