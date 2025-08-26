using CSV
include("../benchmarks/BASBLib_benchmark.jl")
res = benchmark_BASBLib(; do_require_nonstrict_min=true, init_solver="IPOPT", solver="IPOPT");
CSV.write("BASBLib_results.csv", res)

include("../benchmarks/BOLIB_benchmark.jl")
res = benchmark_BOLIB(; do_require_nonstrict_min=true, init_solver="IPOPT", solver="IPOPT");
CSV.write("BOLIB_results.csv", res)

#res_min = run_all_BOLIB_examples(; verbosity=0, max_iter=50, tol=1e-7, conv_dv_len=3, do_check_x_agreem=true, do_force_hp_init=false, max_rand_restart_ct=10, do_require_nonstrict_min=false, init_solver="PATH", solver="PATH");
#df = hcat(DataFrame("SOSC" => fill("", length(res.df[!, 1]))), res.df, DataFrame("SONC" => fill("", length(res.df[!, 1]))), res_min.df[!, 6:12], makeunique=true)
#CSV.write("BOLIB_results.csv", df)
