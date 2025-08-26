if !haskey(ENV, "BOLIB_PATH")
    error("You need to obtain [BOLIB.jl](https://github.com/cinaral/BOLIB.jl) (converted from [BOLIBver2](https://biopt.github.io/bolib/)) and set BOLIB_PATH environment variable to run this benchmark.\n")
end
Pkg.develop(path=ENV["BOLIB_PATH"])
import BOLIB
using BilevelOptSolver
using DataFrames
using Statistics
using Random

"""
Usage:
```
julia> include("benchmarks/BOLIB_benchmark.jl")
julia> df = benchmark_BOLIB(example_ids=1:165);
```
"""
function benchmark_BOLIB(; example_ids=1:length(BOLIB.examples), verbosity=0, tol=1e-7, init_solver="IPOPT", solver="IPOPT", max_iter=50, conv_dv_len=3, do_force_hp_init=false, do_require_nonstrict_min=false, do_check_x_agreem=true, max_rand_restart_ct=50, rng=MersenneTwister(), x_init=[])
    dataframes = []
    success_arr = Bool[]
    elapsed_arr = Float64[]
    prob_count = 0

    for id in example_ids
        name = BOLIB.examples[id]
        print("$id $name\n")

        if "SinhaMaloDeb2014TP9" == name || "SinhaMaloDeb2014TP10" == name
            # these fail to compile for some reason so we skip
            continue
        end

        prob_count += 1
        prob = getfield(Main.BOLIB, Symbol(name))()
        bop, _ = construct_bop(prob.n1, prob.n2, prob.F, prob.G, prob.f, prob.g; verbosity=0)
        if isempty(x_init)
            x_init = gen_x_init(name, rng)
        end

        # dry run for @elapsed...
        is_sol_valid, x, λ, iter_count, status = solve_bop(bop; x_init=prob.xy_init, verbosity=0, tol, init_solver, solver, max_iter=2, conv_dv_len, do_force_hp_init, do_require_nonstrict_min, do_check_x_agreem, max_rand_restart_ct)

        elapsed_time = @elapsed begin
            is_sol_valid, x, λ, iter_count, status = solve_bop(bop; x_init, verbosity, tol, init_solver, solver, max_iter, conv_dv_len, do_force_hp_init, do_require_nonstrict_min, do_check_x_agreem, max_rand_restart_ct)
        end

        Ff = [bop.F(x); bop.f(x)]
        success, rating, x_optimal = rate_BOLIB_result(name, x, Ff, is_sol_valid; tol)

        if is_sol_valid
            info_status = "success"
        else
            info_status = "FAIL"
        end
        info = "$info_status: $rating ($status), $iter_count iters ($(round(elapsed_time, sigdigits=5)) s), x = $(round.(x, sigdigits=5)), Ff = $(round.(Ff, sigdigits=5)) (x* = $(round.(x_optimal, sigdigits=5)), Ff* = $(round.(prob.Ff_optimal, sigdigits=5)))\n"

        if verbosity > 0
            print(info)
        end

        dataframes = [dataframes; DataFrame("name" => name, "n1" => bop.n1, "n2" => bop.n2, "m1" => bop.m1, "m2" => bop.m2, "Iterations" => iter_count, "Status" => status, "Solve time (s)" => round.(elapsed_time, sigdigits=3), "Success" => success, "is_sol_valid" => is_sol_valid, "Rating" => rating, "Ff" => Ref(round.(Ff, sigdigits=3)), "Ff*" => Ref(round.(prob.Ff_optimal, sigdigits=3)), "x" => Ref(round.(x, sigdigits=3)), "x*" => Ref(round.(x_optimal, sigdigits=3)), "x_init" => Ref(round.(x_init, sigdigits=3)))]
        push!(success_arr, success)
        push!(elapsed_arr, elapsed_time)
    end

    success_elapsed_arr = elapsed_arr[success_arr]
    success_count = length(findall(success_arr))

    print("Out of $(prob_count) problems, $(success_count) ($(round((success_count/prob_count*100),sigdigits=3))%) were successful.\n")
    print("Elapsed min-max: $(round(minimum(elapsed_arr),sigdigits=2))-$(round(maximum(elapsed_arr),sigdigits=2)) s, median: $(round(median(elapsed_arr),sigdigits=2)) s, mean: $(round(mean(elapsed_arr),sigdigits=2)) s\n")
    if !isempty(success_elapsed_arr)
        print("Elapsed successful min-max: $(round(minimum(success_elapsed_arr),sigdigits=2))-$(round(maximum(success_elapsed_arr),sigdigits=2)) s, median: $(round(median(success_elapsed_arr),sigdigits=2)) s, mean: $(round(mean(success_elapsed_arr),sigdigits=2))\n")
    end

    vcat(dataframes...)
end

function rate_BOLIB_result(name, x, Ff, is_sol_valid; tol)
    prob = getfield(Main.BOLIB, Symbol(name))()
    Ff_optimal = prob.Ff_optimal[1:2]
    success = false
    rating = ""
    x_optimal = []

    if prob.Ff_optimal[3] == 0 # star means optimal
        rating = "no reference"
    elseif prob.Ff_optimal[3] == 1 # star means optimal
        if isapprox(Ff, Ff_optimal; atol=2 * tol)
            rating = "optimal Ff"
        elseif isapprox(Ff[1], Ff_optimal[1]; atol=2 * tol)
            if Ff[2] < Ff_optimal[2] - tol
                rating = "optimal F and better f"
            else
                rating = "optimal F but worse f"
            end
        elseif Ff[1] < Ff_optimal[1] - tol || (isapprox(Ff[1], Ff_optimal[1]; atol=2 * tol) && Ff[2] < Ff_optimal[2] - tol)
            rating = "better F"
        else
            rating = "SUBOPTIMAL Ff"
        end
    elseif prob.Ff_optimal[3] == 2 # star means best
        if isapprox(Ff, Ff_optimal; atol=2 * tol)
            rating = "best Ff"
        elseif isapprox(Ff[1], Ff_optimal[1]; atol=2 * tol)
            if Ff[2] < Ff_optimal[2] - tol
                rating = "best F and better f"
            else
                rating = "best F but WORSE f"
            end
        elseif Ff[1] < Ff_optimal[1] - tol || (isapprox(Ff[1], Ff_optimal[1]; atol=2 * tol) && Ff[2] < Ff_optimal[2] - tol)
            rating = "better than best Ff"
        else
            rating = "WORSE than best Ff"
        end
    end

    if is_sol_valid
        success = true
    end

    (; success, rating, x_optimal)
end

function gen_x_init(name, rng)
    prob = getfield(Main.BASBLib, Symbol(name))()
    x_init = prob.xy_init

    x_init
end