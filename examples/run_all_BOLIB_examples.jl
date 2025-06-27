# https://biopt.github.io/bolib/
using BilevelOptSolver
using Statistics
include("../src/BOLIB_utils.jl")

res = run_all_BOLIB_examples(; verbosity=1, is_using_PATH=false, is_using_HSL=true);
#res_PATH = run_all_BOLIB_examples(; verbosity=1, is_using_PATH=true, is_using_HSL=true);

success_elapsed_arr = res.elapsed_arr[res.success_arr]

print("Out of $(res.prob_count) problems, $(res.success_count) ($(res.success_count/res.prob_count*100)%) were successful.\n")
print("Out of successful solutions: $(res.optimalish_count) ($(res.optimalish_count/res.success_count*100)%) were optimal or best known, while $(res.suboptimalish_count) ($(res.suboptimalish_count/res.success_count*100)%) were suboptimal or worse than best known.\n")
print("Elapsed min-max: $(minimum(res.elapsed_arr))-$(maximum(res.elapsed_arr)) s, median: $(median(res.elapsed_arr)) s, mean: $(mean(res.elapsed_arr)) s\n")
print("Elapsed successful min-max: $(minimum(success_elapsed_arr))-$(maximum(success_elapsed_arr)) s, median: $(median(success_elapsed_arr)) s, mean: $(mean(success_elapsed_arr)) ")
