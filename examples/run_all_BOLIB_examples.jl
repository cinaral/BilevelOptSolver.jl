# https://biopt.github.io/bolib/
using BilevelOptSolver
using Statistics
using CSV
include("../src/BOLIB_utils.jl")

res = run_all_BOLIB_examples(; verbosity=0, max_iter=200, is_using_HSL=true, tol=1e-6, is_checking_min=true);
#res_no_min = run_all_BOLIB_examples(; verbosity=0, max_iter=200, is_using_HSL=true, tol=1e-6, is_checking_min=false);

@info "is_checking_min TRUE"
success_elapsed_arr = res.elapsed_arr[res.success_arr]
print("Out of $(res.prob_count) problems, $(res.success_count) ($(res.success_count/res.prob_count*100)%) were successful.\n")
print("Out of successful solutions: $(res.optimalish_count) ($(res.optimalish_count/res.success_count*100)%) were optimal or best known, while $(res.suboptimalish_count) ($(res.suboptimalish_count/res.success_count*100)%) were suboptimal or worse than best known.\n")
print("Elapsed min-max: $(minimum(res.elapsed_arr))-$(maximum(res.elapsed_arr)) s, median: $(median(res.elapsed_arr)) s, mean: $(mean(res.elapsed_arr)) s\n")
print("Elapsed successful min-max: $(minimum(success_elapsed_arr))-$(maximum(success_elapsed_arr)) s, median: $(median(success_elapsed_arr)) s, mean: $(mean(success_elapsed_arr)) ")

#@info "is_checking_min FALSE"
#success_elapsed_arr = res.elapsed_arr[res_no_min.success_arr]
#print("Out of $(res_no_min.prob_count) problems, $(res_no_min.success_count) ($(res_no_min.success_count/res_no_min.prob_count*100)%) were successful.\n")
#print("Out of successful solutions: $(res_no_min.optimalish_count) ($(res_no_min.optimalish_count/res_no_min.success_count*100)%) were optimal or best known, while $(res_no_min.suboptimalish_count) ($(res_no_min.suboptimalish_count/res_no_min.success_count*100)%) were suboptimal or worse than best known.\n")
#print("Elapsed min-max: $(minimum(res_no_min.elapsed_arr))-$(maximum(res_no_min.elapsed_arr)) s, median: $(median(res_no_min.elapsed_arr)) s, mean: $(mean(res_no_min.elapsed_arr)) s\n")
#print("Elapsed successful min-max: $(minimum(success_elapsed_arr))-$(maximum(success_elapsed_arr)) s, median: $(median(success_elapsed_arr)) s, mean: $(mean(success_elapsed_arr)) ")

#df = hcat(res.df, DataFrame("non-minimizing" => fill("", length(res.df[!, 1]))), res_no_min.df[!,6:12], makeunique=true)
#CSV.write("BOLIB_results.csv", df)
CSV.write("BOLIB_results.csv", res.df)

