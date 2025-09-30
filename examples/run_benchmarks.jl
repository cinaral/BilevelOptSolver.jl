using CSV

# need
# table order
BASBLib.example_ids=1:length(BASBLib.examples)
BOLIB.example_ids=1:length(BOLIB.examples)

include("../benchmarks/BASBLib_benchmark.jl")
#res = benchmark_BASBLib(; do_require_strict_min=false, init_solver="IPOPT", solver="IPOPT", do_force_dry_run=true);
#CSV.write("BASBLib_ipopt_nonstrict.csv", res)
res = benchmark_BASBLib(; do_require_strict_min=true, init_solver="IPOPT", solver="IPOPT");
res = [res; ]
CSV.write("BASBLib_ipopt.csv", res)
#res = benchmark_BASBLib(; do_require_strict_min=false, init_solver="PATH", solver="PATH", do_force_dry_run=true);
#CSV.write("BASBLib_path_nonstrict.csv", res)
#res = benchmark_BASBLib(; do_require_strict_min=true, init_solver="PATH", solver="PATH");
#CSV.write("BASBLib_path.csv", res)

#include("../benchmarks/BOLIB_benchmark.jl")
#res = benchmark_BOLIB(; do_require_strict_min=false, init_solver="IPOPT", solver="IPOPT", do_force_dry_run=true);
#CSV.write("BOLIB_ipopt_nonstrict.csv", res)
#res = benchmark_BOLIB(; do_require_strict_min=true, init_solver="IPOPT", solver="IPOPT");
#CSV.write("BOLIB_ipopt.csv", res)
#res = benchmark_BOLIB(; do_require_strict_min=false, init_solver="PATH", solver="PATH", do_force_dry_run=true);
#CSV.write("BOLIB_path_nonstrict.csv", res)
#res = benchmark_BOLIB(; do_require_strict_min=true, init_solver="PATH", solver="PATH");
#CSV.write("BOLIB_path.csv", res)
