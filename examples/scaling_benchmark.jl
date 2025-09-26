using BilevelOptSolver
using BenchmarkTools

function RobustPortfolioP1(N; d=2)
    n1::Int64 = N + 1
    n2::Int64 = N

    function F(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        -x[end]
    end

    function G(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        [
            -x[end] + y' * x[1:N];
            x[1:N];
            1.0 - sum(x[1:N]);
            sum(x[1:N]) - 1.0
        ]
    end

    function f(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        y' * x[1:N] - x[end]
    end

    function g(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        i = (1:N)'
        si = ((0.05 / 3 / N) * sqrt.(2 * N * (N + 1) .* i)) .^ d
        yi = 1.15 .+ (0.05 / N) .* i
        [
            sum((abs.(y .- yi)) .^ d ./ si) - 1.5^d;
            y
        ]
    end

    xy_init = ones(n1 + n2)
    Ff_optimal = Float64[-1.15; 0; 2]

    (; n1, n2, F, G, f, g, xy_init, Ff_optimal)
end

function RobustPortfolioP2(N; d=2)
    n1::Int64 = N + 1
    n2::Int64 = N

    function F(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        -x[end]
    end

    function G(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        [
            -x[end] + y' * x[1:N];
            x[1:N];
            1.0 - sum(x[1:N]);
            sum(x[1:N]) - 1.0
        ]
    end

    function f(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        y' * x[1:N] - x[end]
    end

    function g(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        i = (1:N)'
        si = ((0.05 / 3 / N) * sqrt.(2 * N * (N + 1) .* i))
        yi = 1.15 .+ (0.05 / N) .* i
        sx = 1.5 * (1 + sum((x[1:N] .- 1 / N) .^ 2))
        [
            sx^2 - sum((abs.(y .- yi)) .^ 2 ./ si);
            y
        ]
    end

    xy_init = ones(n1 + n2)
    Ff_optimal = Float64[-1.15; 0; 2]

    (; n1, n2, F, G, f, g, xy_init, Ff_optimal)
end

"""
Usage:
```
```
"""
function create_solve(prob; tol=1e-7, verbosity=0, init_solver="PATH", solver="PATH", max_iter=50, conv_dv_len=3, do_force_hp_init=false, do_require_strict_min=true, do_check_x_agreem=true, do_force_toggle=false, max_rand_restart_ct=0)
    bop, _ = construct_bop(prob.n1, prob.n2, prob.F, prob.G, prob.f, prob.g; verbosity=0, np=0)

    function solve(; x_init=prob.xy_init)
        if isempty(x_init)
            x_init = zeros(bop.nx)
        end
        is_sol_valid, x, λ, iter_count, status = solve_bop(bop; x_init, verbosity, tol, init_solver, solver, max_iter, conv_dv_len, do_force_hp_init, do_require_strict_min, do_check_x_agreem, do_force_toggle, max_rand_restart_ct)

        (; is_sol_valid, x, λ, iter_count, status)
    end

    (; solve, bop)
end

# fails to compile
#solve, bop = create_solve(RobustPortfolioP1(50));
#is_sol_valid, x, λ, iter_count, status = solve(x_init=rand(bop.nx));

solve, _ = @btime create_solve(RobustPortfolioP1(10));
bench = @benchmark solve();
include("../benchmarks/RobustPortfolio/10_RobustPortfolioP1.jl")
include("../benchmarks/RobustPortfolio/10_RobustPortfolioP2.jl")
solve_P1, _ = @btime create_solve(RobustPortfolioP1());
solve_P2, _ = @btime create_solve(RobustPortfolioP2());
P1_10_bench = @benchmark solve_P1();
P2_10_bench = @benchmark solve_P2();
include("../benchmarks/RobustPortfolio/20_RobustPortfolioP1.jl")
include("../benchmarks/RobustPortfolio/20_RobustPortfolioP2.jl")
solve_P1, _ = @btime create_solve(RobustPortfolioP1());
solve_P2, _ = @btime create_solve(RobustPortfolioP2());
P1_20_bench = @benchmark solve_P1();
P2_20_bench = @benchmark solve_P2();
include("../benchmarks/RobustPortfolio/30_RobustPortfolioP1.jl")
include("../benchmarks/RobustPortfolio/30_RobustPortfolioP2.jl")
solve_P1, _ = @btime create_solve(RobustPortfolioP1());
solve_P2, _ = @btime create_solve(RobustPortfolioP2());
P1_30_bench = @benchmark solve_P1();
P2_30_bench = @benchmark solve_P2();
include("../benchmarks/RobustPortfolio/40_RobustPortfolioP1.jl")
include("../benchmarks/RobustPortfolio/40_RobustPortfolioP2.jl")
solve_P1, _ = @btime create_solve(RobustPortfolioP1());
solve_P2, _ = @btime create_solve(RobustPortfolioP2());
P1_40_bench = @benchmark solve_P1();
P2_40_bench = @benchmark solve_P2();
include("../benchmarks/RobustPortfolio/50_RobustPortfolioP1.jl")
include("../benchmarks/RobustPortfolio/50_RobustPortfolioP2.jl")
solve_P1, _ = @btime create_solve(RobustPortfolioP1());
solve_P2, _ = @btime create_solve(RobustPortfolioP2());
P1_50_bench = @benchmark solve_P1();
P2_50_bench = @benchmark solve_P2();

include("../benchmarks/RobustPortfolio/60_RobustPortfolioP1.jl")
include("../benchmarks/RobustPortfolio/60_RobustPortfolioP2.jl")
solve_P1, _ = @btime create_solve(RobustPortfolioP1());
solve_P2, _ = @btime create_solve(RobustPortfolioP2());
P1_60_bench = @benchmark solve_P1();
P2_60_bench = @benchmark solve_P2();
include("../benchmarks/RobustPortfolio/70_RobustPortfolioP1.jl")
include("../benchmarks/RobustPortfolio/70_RobustPortfolioP2.jl")
solve_P1, _ = @btime create_solve(RobustPortfolioP1());
solve_P2, _ = @btime create_solve(RobustPortfolioP2());
P1_70_bench = @benchmark solve_P1();
P2_70_bench = @benchmark solve_P2();
include("../benchmarks/RobustPortfolio/80_RobustPortfolioP1.jl")
include("../benchmarks/RobustPortfolio/80_RobustPortfolioP2.jl")
solve_P1, _ = @btime create_solve(RobustPortfolioP1());
solve_P2, _ = @btime create_solve(RobustPortfolioP2());
P1_80_bench = @benchmark solve_P1();
P2_80_bench = @benchmark solve_P2();
include("../benchmarks/RobustPortfolio/90_RobustPortfolioP1.jl")
include("../benchmarks/RobustPortfolio/90_RobustPortfolioP2.jl")
solve_P1, _ = @btime create_solve(RobustPortfolioP1());
solve_P2, _ = @btime create_solve(RobustPortfolioP2());
P1_90_bench = @benchmark solve_P1();
P2_90_bench = @benchmark solve_P2();
include("../benchmarks/RobustPortfolio/100_RobustPortfolioP1.jl")
include("../benchmarks/RobustPortfolio/100_RobustPortfolioP2.jl")
solve_P1, _ = @btime create_solve(RobustPortfolioP1());
solve_P2, _ = @btime create_solve(RobustPortfolioP2());
P1_100_bench = @benchmark solve_P1();
P2_100_bench = @benchmark solve_P2();

P1_benchs = (; N10=P1_10_bench, N20=P1_20_bench, N30=P1_30_bench, N40=P1_40_bench, N50=P1_50_bench, N60=P1_60_bench, N70=P1_70_bench, N80=P1_80_bench, N90=P1_90_bench, N100=P1_100_bench);
P2_benchs = (; N10=P2_10_bench, N20=P2_20_bench, N30=P2_30_bench, N40=P2_40_bench, N50=P2_50_bench, N60=P2_60_bench, N70=P2_70_bench, N80=P2_80_bench, N90=P2_90_bench, N100=P2_100_bench);
#(; is_sol_valid, x, λ, iter_count, status) = @btime solve_P2()

P1_medians = [median(P1_benchs.N10.times); median(P1_benchs.N20.times); median(P1_benchs.N30.times); median(P1_benchs.N40.times); median(P1_benchs.N50.times); median(P1_benchs.N60.times); median(P1_benchs.N70.times); median(P1_benchs.N80.times); median(P1_benchs.N90.times); median(P1_benchs.N100.times)] ./ 1e9

P2_medians = [median(P2_benchs.N10.times); median(P2_benchs.N20.times); median(P2_benchs.N30.times); median(P2_benchs.N40.times); median(P2_benchs.N50.times); median(P2_benchs.N60.times); median(P2_benchs.N70.times); median(P2_benchs.N80.times); median(P2_benchs.N90.times); median(P2_benchs.N100.times)] ./ 1e9

P1_means = [mean(P1_benchs.N10.times); mean(P1_benchs.N20.times); mean(P1_benchs.N30.times); mean(P1_benchs.N40.times); mean(P1_benchs.N50.times); mean(P1_benchs.N60.times); mean(P1_benchs.N70.times); mean(P1_benchs.N80.times); mean(P1_benchs.N90.times); mean(P1_benchs.N100.times)] ./ 1e9

P2_means = [mean(P2_benchs.N10.times); mean(P2_benchs.N20.times); mean(P2_benchs.N30.times); mean(P2_benchs.N40.times); mean(P2_benchs.N50.times); mean(P2_benchs.N60.times); mean(P2_benchs.N70.times); mean(P2_benchs.N80.times); mean(P2_benchs.N90.times); mean(P2_benchs.N100.times)] ./ 1e9

function get_CI(vals; scale=1.0)
    vals = vals .* scale
    CI = 1.96 * std(vals) / sqrt(length(vals))
    CI
end

P1_CI = [get_CI(P1_benchs.N10.times); get_CI(P1_benchs.N20.times); get_CI(P1_benchs.N30.times); get_CI(P1_benchs.N40.times); get_CI(P1_benchs.N50.times); get_CI(P1_benchs.N60.times); get_CI(P1_benchs.N70.times); get_CI(P1_benchs.N80.times); get_CI(P1_benchs.N90.times); get_CI(P1_benchs.N100.times)] ./ 1e9

P2_CI = [get_CI(P2_benchs.N10.times); get_CI(P2_benchs.N20.times); get_CI(P2_benchs.N30.times); get_CI(P2_benchs.N40.times); get_CI(P2_benchs.N50.times); get_CI(P2_benchs.N60.times); get_CI(P2_benchs.N70.times); get_CI(P2_benchs.N80.times); get_CI(P2_benchs.N90.times); get_CI(P2_benchs.N100.times)] ./ 1e9

#P1_means = [0.00347315
#0.0043901
#0.0048312
#0.0062719
#0.0072221
#0.009388
#0.008347
#0.00942525
#0.01040315
#0.0146771
#]
#P2_means = [0.2037161
#0.2309925
#0.67421825
#0.8523594
#1.1406361
#1.5333983
#1.734028
#2.0509055
#2.9242926
#2.8849767
#]

#p = plot(10:10:100, P1_means)
#plot!(p, size=(400, 200), xlabel="N", ylabel="Median " * L"t_c" * " (ms)", yaxis=(formatter = y -> round(y*1e3; sigdigits=4)), legend=false)
#savefig("./output/RobustPorfolioP1_scaling.pdf")

#p = plot(10:10:100, P2_means)
#plot!(p, size=(400, 200), xlabel="N", ylabel="Median " * L"t_c" * " (s)", yaxis=(formatter = y -> round(y; sigdigits=4)), legend=false)
#savefig("./output/RobustPorfolioP2_scaling.pdf")
#plot!(p, 10:10:100, P2_means)

#bop =  bop_P1;
print("")
#bop =  bop_P2 # fails
#include("../src/forrest_solver.jl")
#using .forrest_solver
#OP1 = forrest_solver.OptimizationProblem(bop.n1 + bop.n2, 1:bop.n1, bop.F, bop.G, zeros(bop.m1), fill(Inf, bop.m1))
#OP2 = forrest_solver.OptimizationProblem(bop.n1 + bop.n2, 1:bop.n2, bop.f, bop.g, zeros(bop.m2), fill(Inf, bop.m2))
#bilevel = [OP1; OP2]
#out_bilevel = @btime forrest_solver.solve(bilevel)
#x_forrest = out_bilevel[1:bop.n1+bop.n2]
#@info x_forrest
#@info bop.F(x_forrest)

#Ff = [bop.F(x); bop.f(x)]
#rating = rate_BASBLib_result(name, x, Ff; tol=1e-3)
#if is_sol_valid
#    info_status = "success"
#else
#    info_status = "FAIL"
#end
#info = "$info_status: $rating ($status), $iter_count iters ($(round(elapsed_time, sigdigits=5)) s), x = $(round.(x, sigdigits=5)), Ff = $(round.(Ff, sigdigits=5)) (x* = $(round.(prob.xy_optimal, sigdigits=5)), Ff* = $(round.(prob.Ff_optimal, sigdigits=5)))"

#if verbosity > 0
#    print(info)
#end

#(; info, Ff, is_sol_valid, x, λ, iter_count, status, elapsed_time, bop, syms, solve)