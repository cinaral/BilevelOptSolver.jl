using CSV
#include("../benchmarks/BASBLib_benchmark.jl")
#res = benchmark_BASBLib(; do_require_strict_min=false, init_solver="IPOPT", solver="IPOPT", do_force_dry_run=true);
#CSV.write("BASBLib_results_nonstrict.csv", res)
#res = benchmark_BASBLib(; do_require_strict_min=true, init_solver="IPOPT", solver="IPOPT");
#CSV.write("BASBLib_results.csv", res)
#res = benchmark_BASBLib(; do_require_strict_min=false, init_solver="PATH", solver="PATH", do_force_dry_run=true);
#CSV.write("BASBLib_results_path_nonstrict.csv", res)
#res = benchmark_BASBLib(; do_require_strict_min=true, init_solver="PATH", solver="PATH");
#CSV.write("BASBLib_results_path.csv", res)

include("../benchmarks/BOLIB_benchmark.jl")
#res = benchmark_BOLIB(; do_require_strict_min=false, init_solver="IPOPT", solver="IPOPT", do_force_dry_run=true);
#CSV.write("BOLIB_results_nonstrict.csv", res)
res = benchmark_BOLIB(; do_require_strict_min=true, init_solver="IPOPT", solver="IPOPT");
CSV.write("BOLIB_results.csv", res)
#res = benchmark_BOLIB(; do_require_strict_min=false, init_solver="PATH", solver="PATH", do_force_dry_run=true);
#CSV.write("BOLIB_results_path_nonstrict.csv", res)
#res = benchmark_BOLIB(; do_require_strict_min=true, init_solver="PATH", solver="PATH");
#CSV.write("BOLIB_results_path.csv", res)
