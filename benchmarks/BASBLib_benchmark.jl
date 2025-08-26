
if !haskey(ENV, "BASBLib_PATH")
    error("You need to obtain [BASBLib.jl](https://github.com/cinaral/BASBLib.jl) (converted from [BASBLibv2.2](https://github.com/basblsolver/BASBLib)) and set BASBLib_PATH environment variable to run this benchmark.\n")
end
Pkg.develop(path=ENV["BASBLib_PATH"])
import BASBLib
using BilevelOptSolver
using DataFrames
using Statistics
using Random

"""
Usage:
```
julia> include("benchmarks/BASBLib_benchmark.jl")
julia> df = benchmark_BASBLib(example_ids=1:82);
```
"""
function benchmark_BASBLib(; example_ids=1:length(BASBLib.examples), verbosity=0, tol=1e-7, init_solver="IPOPT", solver="IPOPT", max_iter=50, conv_dv_len=3, do_force_hp_init=false, do_require_nonstrict_min=false, do_check_x_agreem=true, max_rand_restart_ct=50, rng=MersenneTwister(), x_init=[], max_retry_ct=1)
    dataframes = []
    success_arr = Bool[]
    elapsed_arr = Float64[]
    prob_count = 0

    for id in example_ids
        name = BASBLib.examples[id]
        print("$id $name\n")

        prob_count += 1
        prob = getfield(Main.BASBLib, Symbol(name))()
        bop, _ = construct_bop(prob.n1, prob.n2, prob.F, prob.G, prob.f, prob.g; verbosity=0)
        if isempty(x_init)
            x_init = gen_x_init(name, rng)
        end

        # dry run for @elapsed and can be used to check for init issues
        is_sol_valid = false
        retry_ct = 0
        while !is_sol_valid && retry_ct < max_retry_ct
            @info x_init
            is_sol_valid, x, λ, iter_count, status = solve_bop(bop; x_init, verbosity=0, tol, init_solver, solver, max_iter, conv_dv_len, do_force_hp_init, do_require_nonstrict_min, do_check_x_agreem, max_rand_restart_ct)
            if !is_sol_valid
                @info "poop"
                x_init = gen_x_init(name, rng)
            end
            retry_ct += 1
        end

        elapsed_time = @elapsed begin
            is_sol_valid, x, λ, iter_count, status = solve_bop(bop; x_init, verbosity, tol, init_solver, solver, max_iter, conv_dv_len, do_force_hp_init, do_require_nonstrict_min, do_check_x_agreem, max_rand_restart_ct)
        end

        Ff = [bop.F(x); bop.f(x)]
        success, rating = rate_BASBLib_result(name, x, Ff, is_sol_valid; tol)
        if is_sol_valid
            info_status = "success"
        else
            info_status = "FAIL"
        end
        info = "$info_status: $rating ($status), $iter_count iters ($(round(elapsed_time, sigdigits=5)) s), x = $(round.(x, sigdigits=5)), Ff = $(round.(Ff, sigdigits=5)) (x* = $(round.(prob.xy_optimal, sigdigits=5)), Ff* = $(round.(prob.Ff_optimal, sigdigits=5)))\n"

        if verbosity > 0
            print(info)
        end

        dataframes = [dataframes; DataFrame("name" => name, "n1" => bop.n1, "n2" => bop.n2, "m1" => bop.m1, "m2" => bop.m2, "Iterations" => iter_count, "Status" => status, "Solve time (s)" => round.(elapsed_time, sigdigits=3), "Success" => success, "is_sol_valid" => is_sol_valid, "Rating" => rating, "Ff" => Ref(round.(Ff, sigdigits=3)), "Ff*" => Ref(round.(prob.Ff_optimal, sigdigits=3)), "x" => Ref(round.(x, sigdigits=3)), "x*" => Ref(round.(prob.xy_optimal, sigdigits=3)), "x_init" => Ref(round.(x_init, sigdigits=3)))]
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

function rate_BASBLib_result(name, x, Ff, is_sol_valid; tol=1e-7)
    prob = getfield(Main.BASBLib, Symbol(name))()
    success = false
    rating = ""

    if isempty(prob.Ff_optimal)
        rating = "no reference"
    else
        is_cost_optimal = isapprox(Ff, prob.Ff_optimal; atol=tol)
        is_xy_optimal = isapprox(x, prob.xy_optimal; atol=tol)

        if (is_cost_optimal && is_xy_optimal)
            rating = "optimal"
        else
            rating = "NOT optimal"
        end
    end

    if name == "mb_2007_02" # no solution
        if !is_sol_valid
            success = true
        end
    else
        if rating == "optimal"
            success = true
        end
    end

    (; success, rating)
end

function gen_x_init(name, rng)
    prob = getfield(Main.BASBLib, Symbol(name))()
    x_init = zeros(prob.n1 + prob.n2)

    if name == "mb_2007_02"
        x_init[1] = -rand(rng)
    elseif name == "lh_1994_01"
        x_init = rand(rng, Float64, 2) .* 10.0
    end
    x_init
end